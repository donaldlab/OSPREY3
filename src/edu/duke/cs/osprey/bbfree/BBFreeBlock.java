/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.bbfree;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DOFBlock;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.deeper.GenChi1Calc;
import edu.duke.cs.osprey.dof.deeper.SidechainIdealizer;
import edu.duke.cs.osprey.ematrix.epic.SeriesFitter;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.RigidBodyMotion;
import edu.duke.cs.osprey.tools.RotationMatrix;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

/**
 *
 * This is a set of consecutive residues whose backbones can move freely,
 * subject to distance, angle, and omega constraints
 * The freely moving zone is anchored fix CA's at either end
 * 
 * @author mhall44
 */
public class BBFreeBlock implements Serializable, DOFBlock {
    
    List<Residue> residues;//the residues moving freely, in order 
    //including the two end ones (only carboxyl or amine moving freely)
    
    ArrayList<BBFreeDOF> freeDOFs = new ArrayList<>();
    
    double[][] freeDOFVoxel;//voxel of allowed values for free DOFs.  We'll minimize over this.
    //Indices: 0/1 (for min/max), DOF # 
    
    double[] curFreeDOFVals;//current values, so we can update one at a time

    double[] fullDOFCenter;//values of full DOFs at center of voxel
    DoubleMatrix1D freeDOFCenter;
    
    double[][] fullDOFPolys;//Polynomials represented in double[] format (see SeriesFitter)
    
    
    PepPlaneLinModel[] pepPlanes;//linear models of BB atoms in peptide planes
    //as a function of N and CA coordinates
    
    
    DoubleMatrix2D freeDOFMatrix;//matrix whose rows are the coefficients defining free DOFs
    
    static int polyOrder = 2;
    
    
    ArrayList<ArrayList<JacDerivEntry>> jacDerivs;
    
    
    public BBFreeBlock(){
        
    }
    
    public BBFreeBlock(List<Residue> residues){
        // Calculate the free DOFs, indicate appropriate voxels for them
        // and fit polynomials to express the full set of polynomials as a function of them
        
        this.residues = residues;
        int numRes = residues.size();
        int numFreeDOFs = 2*numRes-6;
        
        pepPlanes = new PepPlaneLinModel[numRes-1];
        for(int planeNum=0; planeNum<numRes-1; planeNum++)
            pepPlanes[planeNum] = new PepPlaneLinModel(residues.get(planeNum),residues.get(planeNum+1));
        
        //generating free DOFs based on current geometry 
        DoubleMatrix2D constrJac = getConstrJac();//Evaluate at current CA, N coords
        DoubleMatrix2D freeDOFCoeffs = getOrthogVectors(constrJac);
        
        freeDOFMatrix = Algebra.DEFAULT.transpose(freeDOFCoeffs);
        
        for(int freeDOF=0; freeDOF<numFreeDOFs; freeDOF++){
            freeDOFs.add( new BBFreeDOF(freeDOFCoeffs.viewColumn(freeDOF), this, freeDOFs.size()) );
        }
        
        fullDOFCenter = new double[6*numRes-9];
        
        //Pick voxels?  START W/ UNIT CUBE in middle fading quadratically
        //The function 1 - 4(x-0.5)^2 goes from 0 to 1 to 0 on the unit interval
        //the "anchoring" fully fixed res are effectively residues # -1 and numRes
        for(int resNum=1; resNum<numRes; resNum++){
            
            Residue curRes = residues.get(resNum);
            //set init coords
            System.arraycopy(curRes.getCoordsByAtomName("N"), 0, fullDOFCenter, 6*(resNum-1), 3);
            if(resNum<numRes-1)//CA free too
                System.arraycopy(curRes.getCoordsByAtomName("CA"), 0, fullDOFCenter, 6*resNum-3, 3);
        }        

        freeDOFCenter = DoubleFactory1D.dense.make(numFreeDOFs);
        for(int f=0; f<numFreeDOFs; f++)
            freeDOFCenter.set( f, freeDOFs.get(f).evalAtFullDOFs(DoubleFactory1D.dense.make(fullDOFCenter)) );
        
        DoubleMatrix2D J = DoubleFactory2D.dense.compose(
                new DoubleMatrix2D[][] { new DoubleMatrix2D[] {Algebra.DEFAULT.transpose(freeDOFCoeffs)},
                    new DoubleMatrix2D[] {constrJac}
                } );//Jacobian of free DOFs and then constr vars
        makeTaylorSeries(J);
        
        setVoxelBySeries();//set voxel, adjust size if needed to ensure series validity
        //to specified residual in constraints at indicated fullDOF vals (maybe 0.01 distance err?)
        
        curFreeDOFVals = new double[numFreeDOFs];//freeDOFCenter.copy().toArray();//may move away from center, so copy
        //treating these as relative!  for voxel and setDOFs purposes
    }
    
    private void setVoxelBySeries(){
        int numRes = residues.size();
        int numFreeDOFs = 2*numRes-6;
        int numFullDOFs = 6*numRes-9;
        
        freeDOFVoxel = new double[2][numFreeDOFs];
        Arrays.fill(freeDOFVoxel[0],-1);
        Arrays.fill(freeDOFVoxel[1],1);
        //THESE ARE RELATIVE TO CENTER
        
        double sizeStepFactor = 1.3;//how much we'll try decreasing it at a time, if not good enough
        
        double targetResid = 1e-4;//max residual allowed in bond lengths or dot products
        VoxelSeriesChecker vsc = new VoxelSeriesChecker(residues,numFreeDOFs,numFullDOFs,
            fullDOFPolys, pepPlanes, freeDOFMatrix, freeDOFCenter);
        
        
        double resid = vsc.getConstraintsResid(freeDOFVoxel);
        System.out.println("Init voxel resid: "+resid);
        
        while (resid>targetResid) {
            
            for(int a=0; a<2; a++){
                for(int f=0; f<numFreeDOFs; f++)
                    freeDOFVoxel[a][f] /= sizeStepFactor;
            }

            resid = vsc.getConstraintsResid(freeDOFVoxel);
            
            System.out.println("Voxel reduced resid: "+resid);
        }
    }
    
    
       
    //Cached stuff for derivative calcs
    //i, v, w will be full-DOF ("x") indices
    //u, j, k, l, m are free-DOF indices (we ultimately differentiate wrt these)
    //If > numFreeDOFs is used as a constr-DOF index
    double[][] deriv1;//Indices: i, j.  dx_i/dy_j
    double[][][] wderiv2;//Indices: i, j, w.  d^2 x_i / dy_j dx_w
    double[][][] deriv2;//d^x x_i / dy_j dy_k
    double[][][][] wderiv3;//d^3 x_i / dy_j dx_w dy_l
    double[][][][] deriv3;//dx^4 x_i / dy_j dy_k dy_l
    double[][][][][] wderiv4;//d^4 x_i / dy_j dx_w dy_l dy_m
    double[][][][][] deriv4;//d^4 x_i / dy_j dy_k dy_l dy_m
    
    
    private void allocateDerivs(int numFullDOFs, int numFreeDOFs){
        deriv1 = new double[numFullDOFs][numFullDOFs];
        wderiv2 = new double[numFullDOFs][numFullDOFs][numFullDOFs];
        deriv2 = new double[numFullDOFs][numFullDOFs][numFreeDOFs];
        wderiv3 = new double[numFullDOFs][numFullDOFs][numFullDOFs][numFreeDOFs];
        deriv3 = new double[numFullDOFs][numFullDOFs][numFreeDOFs][numFreeDOFs];
        wderiv4 = new double[numFullDOFs][numFullDOFs][numFullDOFs][numFreeDOFs][numFreeDOFs];
        deriv4 = new double[numFullDOFs][numFullDOFs][numFreeDOFs][numFreeDOFs][numFreeDOFs];
    }
    
    private void calcDerivs(int numFullDOFs, int numFreeDOFs, DoubleMatrix2D Jinv){
        //Calculate and cache derivatives for use in Taylor series
        
        allocateDerivs(numFullDOFs, numFreeDOFs);
        
        //1st order
        for(int i=0; i<numFullDOFs; i++){
            for(int j=0; j<numFullDOFs; j++)
                deriv1[i][j] = Jinv.get(i, j);
        }
        
        //2nd order
        for(int i=0; i<numFullDOFs; i++){
            for(int j=0; j<numFullDOFs; j++){
                for(int w=0; w<numFullDOFs; w++){
                    wderiv2[i][j][w] = 0;
                    for(JacDerivEntry jde : jacDerivs.get(w))
                        wderiv2[i][j][w] -= jde.val * deriv1[i][jde.u]  * deriv1[jde.v][j];
                }
                
                for(int k=0; k<numFreeDOFs; k++){
                    deriv2[i][j][k] = 0;
                    for(int w=0; w<numFullDOFs; w++)
                        deriv2[i][j][k] += wderiv2[i][j][w] * deriv1[w][k];
                }
            }
        }
        
        //3rd order
        if(polyOrder>=3){
            for(int i=0; i<numFullDOFs; i++){
                for(int j=0; j<numFullDOFs; j++){
                    
                    for(int l=0; l<numFreeDOFs; l++){
                    
                        for(int w=0; w<numFullDOFs; w++){
                            wderiv3[i][j][w][l] = 0;
                            for(JacDerivEntry jde : jacDerivs.get(w)){
                                wderiv3[i][j][w][l] -= jde.val * deriv2[i][jde.u][l]  * deriv1[jde.v][j];
                                wderiv3[i][j][w][l] -= jde.val * deriv1[i][jde.u]  * deriv2[jde.v][j][l];
                            }
                        }

                        for(int k=0; k<numFreeDOFs; k++){
                            deriv3[i][j][k][l] = 0;
                            for(int w=0; w<numFullDOFs; w++){
                                deriv3[i][j][k][l] += wderiv3[i][j][w][l] * deriv1[w][k];
                                deriv3[i][j][k][l] += wderiv2[i][j][w] * deriv2[w][k][l];
                            }
                        }
                    }
                }
            }
        }
        
        //4th order
        if(polyOrder>=4){
            for(int i=0; i<numFullDOFs; i++){
                for(int j=0; j<numFullDOFs; j++){
                    
                    for(int l=0; l<numFreeDOFs; l++){
                        for(int m=0; m<numFreeDOFs; m++){
                    
                            for(int w=0; w<numFullDOFs; w++){
                                wderiv4[i][j][w][l][m] = 0;
                                for(JacDerivEntry jde : jacDerivs.get(w)){
                                    wderiv4[i][j][w][l][m] -= jde.val * deriv3[i][jde.u][l][m]  * deriv1[jde.v][j];
                                    wderiv4[i][j][w][l][m] -= jde.val * deriv2[i][jde.u][l]  * deriv2[jde.v][j][m];
                                    wderiv4[i][j][w][l][m] -= jde.val * deriv2[i][jde.u][m]  * deriv2[jde.v][j][l];
                                    wderiv4[i][j][w][l][m] -= jde.val * deriv1[i][jde.u]  * deriv3[jde.v][j][l][m];
                                }
                            }

                            for(int k=0; k<numFreeDOFs; k++){
                                deriv4[i][j][k][l][m] = 0;
                                for(int w=0; w<numFullDOFs; w++){
                                    deriv4[i][j][k][l][m] += wderiv4[i][j][w][l][m] * deriv1[w][k];
                                    deriv4[i][j][k][l][m] += wderiv3[i][j][w][l] * deriv2[w][k][m];
                                    deriv4[i][j][k][l][m] += wderiv3[i][j][w][m] * deriv2[w][k][l];
                                    deriv4[i][j][k][l][m] += wderiv2[i][j][w] * deriv3[w][k][l][m];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    
    
    private void makeTaylorSeries(DoubleMatrix2D J){
        //Build fullDOFPolys as Taylor series, centered at orig conf
        //J is d(freeDOFs,constrVars)/d(fullDOFs)
        
        
        int numFullDOFs = 6*residues.size() - 9;
        int numFreeDOFs = 2*residues.size() - 6;
        int numParams = SeriesFitter.getNumParams(numFreeDOFs, true, polyOrder);
        
        DoubleMatrix2D Jinv = Algebra.DEFAULT.inverse(J);//d(fullDOFs)/d(freeDOFs,constrVars)
        
        calcDerivs(numFullDOFs, numFreeDOFs, Jinv);
        
        fullDOFPolys = new double[numFullDOFs][];
        
        for(int fullDOF=0; fullDOF<numFullDOFs; fullDOF++){
            
            double[] series = new double[numParams];
            series[0] = fullDOFCenter[fullDOF];
            
            for(int freeDOF=0; freeDOF<numFreeDOFs; freeDOF++){
                series[freeDOF+1] = deriv1[fullDOF][freeDOF];
            }
            
            int coeffCount = numFreeDOFs+1;
            
            if(polyOrder>=2){

                for(int dof=0; dof<numFreeDOFs; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){

                        series[coeffCount] = deriv2[fullDOF][dof][dof2] / multFacProd(dof,dof2);
                        coeffCount++;
                    }
                }
            }
            
            if(polyOrder>=3){
                for(int dof=0; dof<numFreeDOFs; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                        for(int dof3=0; dof3<=dof2; dof3++){
                            
                            series[coeffCount] = deriv3[fullDOF][dof][dof2][dof3] / multFacProd(dof,dof2,dof3);
                            coeffCount++;
                        }
                    }
                }
            }
            
            if(polyOrder>=4){
                for(int dof=0; dof<numFreeDOFs; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                        for(int dof3=0; dof3<=dof2; dof3++){
                            for(int dof4=0; dof4<=dof3; dof4++){
                            
                                series[coeffCount] = deriv4[fullDOF][dof][dof2][dof3][dof4] / multFacProd(dof,dof2,dof3,dof4);
                                coeffCount++;
                            }
                        }
                    }
                }
            }
            
                    
            fullDOFPolys[fullDOF] = series;
        }
    }
    
    
    private static int multFacProd(int... a){
        //product of factorials of multiplicites of integers
        //argument will be in nonascending order
        //this is the denominator we need for Taylor series terms (enumerated w/o repeats)
        int curFactor = 1;
        int ans = 1;
        
        for(int b=0; b<a.length-1; b++){
            if(a[b]==a[b+1])
                curFactor++;
            else
                curFactor = 1;
            
            ans *= curFactor;
        }
        
        return ans;
    }

    @Override
    public DOFBlock copyForNewMolecule(Molecule mol, LinkedHashMap<DegreeOfFreedom, DegreeOfFreedom> copiedDOFMap) {
        BBFreeBlock copiedBlock = new BBFreeBlock();
        
        copiedBlock.residues = new ArrayList<>();
        for(Residue res : residues){
            Residue newMolRes = mol.getResByPDBResNumber(res.getPDBResNumber());
            copiedBlock.residues.add(newMolRes);
        }
        
        copiedBlock.freeDOFs = new ArrayList<>();
        for(BBFreeDOF dof : freeDOFs){
            BBFreeDOF copiedDOF = new BBFreeDOF(dof.coeffs, copiedBlock, dof.indexInBlock);
            copiedDOFMap.put(dof, copiedDOF);
            copiedBlock.freeDOFs.add(copiedDOF);
        }
        
        //shallow copy is OK for read-only fields
        copiedBlock.freeDOFVoxel = freeDOFVoxel;
        copiedBlock.curFreeDOFVals = curFreeDOFVals.clone();
        copiedBlock.fullDOFCenter = fullDOFCenter;
        copiedBlock.freeDOFCenter = freeDOFCenter;
        copiedBlock.fullDOFPolys = fullDOFPolys;
        copiedBlock.pepPlanes = pepPlanes;
        copiedBlock.freeDOFMatrix = freeDOFMatrix;
        copiedBlock.jacDerivs = jacDerivs;
        //no need to copy all the deriv1, deriv2, etc. since they're just used to construct the block
        
        return copiedBlock;
    }
    
    
    
    private class JacDerivEntry implements Serializable {
        //The derivative of the Jacobian with respect to the full DOFs is sparse
        //Each of these entries, when placed in jacDerivs[w],
        //indicates that d^2 (free DOF u)/ d(full DOF v) d(full DOF w) = val
        int u, v;
        double val;

        public JacDerivEntry(int u, int v, double val) {
            this.u = u;
            this.v = v;
            this.val = val;
        }
    }
    
    
    private static double invMatrixDeriv(DoubleMatrix2D Minv, int i, int j, int u, int v){
        //d(M^-1)_ij / dM_uv
        return - Minv.get(i, u) * Minv.get(v, j);
    }
    
    private DoubleMatrix2D getOrthogVectors(DoubleMatrix2D M){
        //Get a matrix whose columns are orthogonal to the row space of M
        //which expected to be nonsingular
        DoubleMatrix2D Maug = DoubleFactory2D.dense.make(M.columns(), M.columns());
        Maug.viewPart(0, 0, M.rows(), M.columns()).assign(M);
        
        SingularValueDecomposition svd = new SingularValueDecomposition(Maug);
        
        int numOrthVecs = M.columns() - M.rows();
        if(svd.rank() != M.rows()){
            throw new RuntimeException("ERROR: Singularity in constr jac.  Rank: "+svd.rank());
        }
        
        DoubleMatrix2D orthVecs = svd.getV().viewPart(0, M.rows(), M.columns(), numOrthVecs);
        
        //DEBUG!!!  Should be 0 and identity respecitvely
        /*DoubleMatrix2D check = Algebra.DEFAULT.mult(M, orthVecs);
        DoubleMatrix2D checkOrth = Algebra.DEFAULT.mult( Algebra.DEFAULT.transpose(orthVecs), orthVecs);
        */
        
        
        return orthVecs;
    }
    
    
    private DoubleMatrix2D getConstrJac(){
        //Jacobian constraints, evaluated at current (original) coordinates
        
        int numRes = residues.size();
        int numFullDOFs = 6*numRes-9;
        int numFreeDOFs = 2*numRes-6;
        
        jacDerivs = new ArrayList<ArrayList<JacDerivEntry>>();
        for(int f=0; f<numFullDOFs; f++)
            jacDerivs.add(new ArrayList<JacDerivEntry>());
        
        DoubleMatrix2D ans = DoubleFactory2D.dense.make(numFullDOFs-numFreeDOFs,numFullDOFs);
        
        
        int constrCount = 0;
        
        for(int resNum=0; resNum<numRes; resNum++){
            
            Residue curRes = residues.get(resNum);
            
            double NCoord[] = curRes.getCoordsByAtomName("N");
            double CACoord[] = curRes.getCoordsByAtomName("CA");
            
            if(resNum>0){//these constrs act only on fixed coords for resNum=0
                
                Residue prevRes = residues.get(resNum-1);
                double prevCACoord[] = prevRes.getCoordsByAtomName("CA");
                
                //N to CA constr for resNum
                makeDistConstrJac(ans,constrCount,NCoord,CACoord,2*(resNum-1),2*resNum-1);
                constrCount++;

                //CA to CA constr
                makeDistConstrJac(ans,constrCount,prevCACoord,CACoord,2*resNum-3,2*resNum-1);
                constrCount++;

                //last CA to N constr
                makeDistConstrJac(ans,constrCount,prevCACoord,NCoord,2*resNum-3,2*(resNum-1));
                constrCount++;
            }
            
            //OK and finally N-CA-C' angle constr
            if(resNum==numRes-1){//C' fixed: special form of constraint
                double CCoord[] = curRes.getCoordsByAtomName("C");
                makeLastNCACCConstrJac( ans, constrCount, CACoord, CCoord, 2*(resNum-1) );
            }
            else {
                Residue nextRes = residues.get(resNum+1);
                double NCoordNext[] = nextRes.getCoordsByAtomName("N");
                double CACoordNext[] = nextRes.getCoordsByAtomName("CA");

                makeNCACCConstrJac(ans, constrCount, NCoord, CACoord, 
                    NCoordNext, CACoordNext, 2*(resNum-1), resNum);
            }
            
            constrCount++;
        }
                
        return ans;
    }
    
    
    private void recordJacDeriv(int constrNum, int fullDOF1, int fullDOF2, double val){
        int numFreeDOFs = 2*residues.size()-6;
        jacDerivs.get(fullDOF1).add( new JacDerivEntry(constrNum+numFreeDOFs,fullDOF2,val) );
    }
    
    
    private void makeDistConstrJac( DoubleMatrix2D constrJac, int constrNum,
            double coord1[], double coord2[],
            int atomNum1, int atomNum2 ) {
        //Enter the Jacobian of a distance constraint between the atoms as #constrNum
        //in constrJac.  Atom numbers (among free N, CA) and (initial) coordinates given
        //-1 for atomNum means fixed atom, so not part of gradient
        
        for(int dim=0; dim<3; dim++){
            if(atomNum1>=0 && 3*atomNum1<constrJac.columns()){
                constrJac.set(constrNum, 3*atomNum1+dim, 2*coord1[dim]-2*coord2[dim]);
                
                recordJacDeriv(constrNum, 3*atomNum1+dim, 3*atomNum1+dim, 2);
                if(atomNum2>=0 && 3*atomNum2<constrJac.columns())//atomNum2 free to move
                    recordJacDeriv(constrNum, 3*atomNum1+dim, 3*atomNum2+dim, -2);
            }
            
            if(atomNum2>=0 && 3*atomNum2<constrJac.columns()){
                constrJac.set(constrNum, 3*atomNum2+dim, 2*coord2[dim]-2*coord1[dim]);
                
                recordJacDeriv(constrNum, 3*atomNum2+dim, 3*atomNum2+dim, 2);
                if(atomNum1>=0 && 3*atomNum1<constrJac.columns())//atomNum1 free to move
                    recordJacDeriv(constrNum, 3*atomNum2+dim, 3*atomNum1+dim, -2);
            }
        }
    }
    
    
    private void makeNCACCConstrJac( DoubleMatrix2D constrJac, int constrNum,
            double NCoord[], double CACoord[], double NCoordNext[], double CACoordNext[],
            int atomNumN, int pepPlaneNum ){ 
        //constraint the dot product (C'-CA) dot (N-CA)
        //Uses N, CA coord variables for this residue and the next
        //(since C' coords are a linear function of this CA and next N, CA coords)
        //(For this purpose we constrain not the actual C' but its projection into the 
        //CA-N-C plane)
        
        //we'll need expansion coefficients for C'...
        double[] expansionCoeffs = pepPlanes[pepPlaneNum].getProjCAtomCoeffs();
        double a = expansionCoeffs[0];
        double b = expansionCoeffs[1];
        double c = expansionCoeffs[2];
        
        for(int dim=0; dim<3; dim++){
            if(atomNumN>=0){//not the first residue, which has fixed N & CA
                double NGrad = (a-1)*CACoord[dim] + b*NCoordNext[dim] + c*CACoordNext[dim];
                double CAGrad  = (a-1)*NCoord[dim] + 2*(1-a)*CACoord[dim] - b*NCoordNext[dim] - c*CACoordNext[dim];
                constrJac.set(constrNum, 3*atomNumN+dim, NGrad);
                constrJac.set(constrNum, 3*atomNumN+3+dim, CAGrad);
                
                //jac derivs for N
                recordJacDeriv(constrNum, 3*atomNumN+dim, 3*atomNumN+3+dim, a-1);
                recordJacDeriv(constrNum, 3*atomNumN+dim, 3*atomNumN+6+dim, b);
                if(3*atomNumN+9<constrJac.columns())
                    recordJacDeriv(constrNum, 3*atomNumN+dim, 3*atomNumN+9+dim, c);
                
                //and for C
                recordJacDeriv(constrNum, 3*atomNumN+3+dim, 3*atomNumN+dim, a-1);
                recordJacDeriv(constrNum, 3*atomNumN+3+dim, 3*atomNumN+3+dim, 2*(1-a));
                recordJacDeriv(constrNum, 3*atomNumN+3+dim, 3*atomNumN+6+dim, -b);
                if(3*atomNumN+9<constrJac.columns())
                    recordJacDeriv(constrNum, 3*atomNumN+3+dim, 3*atomNumN+9+dim, -c);
            }
            
            constrJac.set(constrNum, 3*atomNumN+6+dim, b*NCoord[dim]-b*CACoord[dim]);
            
            if(atomNumN>=0){
                recordJacDeriv(constrNum, 3*atomNumN+6+dim, 3*atomNumN+dim, b);
                recordJacDeriv(constrNum, 3*atomNumN+6+dim, 3*atomNumN+3+dim, -b);
            }
            
            
            if(3*atomNumN+9<constrJac.columns()){//i.e., not second-to-last res (which has fixed CANext)
                constrJac.set(constrNum, 3*atomNumN+9+dim, c*NCoord[dim]-c*CACoord[dim]);
                
                if(atomNumN>=0){
                    recordJacDeriv(constrNum, 3*atomNumN+9+dim, 3*atomNumN+dim, c);
                    recordJacDeriv(constrNum, 3*atomNumN+9+dim, 3*atomNumN+3+dim, -c);
                }
            }
        }
    }
    
    
    private void makeLastNCACCConstrJac( DoubleMatrix2D constrJac, int constrNum,
            double CACoord[], double CCoord[], int atomNumN ){ 
        //Last constraint is simpler, because CA and C' are both fixed
        
        for(int dim=0; dim<3; dim++)
            constrJac.set(constrNum, 3*atomNumN+dim, CCoord[dim] - CACoord[dim]);
        
        //no jac derivs because this constraint is linear
    }
    
    
    public void setDOFs(DoubleMatrix1D x){
        //x: free DOFs (relative to center, so can eval polys directly)
        
        int numRes = residues.size();
        
        int numDOFsFull = fullDOFPolys.length;
        double fullDOFVals[] = new double[fullDOFPolys.length];
        
        for(int fullDOF=0; fullDOF<numDOFsFull; fullDOF++)
            fullDOFVals[fullDOF] = SeriesFitter.evalSeries(fullDOFPolys[fullDOF], x, x.size(), true, polyOrder);
                    
        //record current information needed for placement of sidechain
        double genChi1[] = new double[numRes];
        double SCTranslations[][] = new double[numRes][3];
        
        //now set each residue in the right place
        
        //start with the CA's and N's
        for(int resNum=0; resNum<numRes; resNum++){
            
            Residue curRes = residues.get(resNum);
            
            //the sidechain will be translated based on CA motion
            if(resNum>0 && resNum<numRes-1){//CA moves, so translate sidechain with it
                double curCACoords[] = curRes.getCoordsByAtomName("CA");
                for(int dim=0; dim<3; dim++)
                    SCTranslations[resNum][dim] = fullDOFVals[6*(resNum-1)+3+dim] - curCACoords[dim];
            }
            
            genChi1[resNum] = GenChi1Calc.getGenChi1(curRes);
            
            //OK now place the backbone...
            if(resNum>0){//N or CA moves
                int CAIndex = curRes.getAtomIndexByName("CA");
                int NIndex = curRes.getAtomIndexByName("N");

                for(int dim=0; dim<3; dim++){
                    if(resNum<numRes-1)//CA moves
                        curRes.coords[3*CAIndex+dim] = fullDOFVals[6*(resNum-1)+3+dim];
                    
                    curRes.coords[3*NIndex+dim] = fullDOFVals[6*(resNum-1)+dim];
                }
            }
        }
        
        //OK now that the CA's and N's are in place, we can finish each peptide plane
        for(int pepPlaneNum=0; pepPlaneNum<numRes-1; pepPlaneNum++){
            //set each atom based on CA and N (linear relationship)
            Residue res1 = residues.get(pepPlaneNum);//residue at start of peptide plane...
            Residue res2 = residues.get(pepPlaneNum+1);//...and at end
            
            double CA1[] = res1.getCoordsByAtomName("CA");
            double N[] = res2.getCoordsByAtomName("N");
            double CA2[] = res2.getCoordsByAtomName("CA");
                        
            res1.setCoordsByAtomName("C", pepPlanes[pepPlaneNum].calcCCoords(CA1,N,CA2,false));
            res1.setCoordsByAtomName("O", pepPlanes[pepPlaneNum].calcOCoords(CA1,N,CA2));
            res2.setCoordsByAtomName("H", pepPlanes[pepPlaneNum].calcHCoords(CA1,N,CA2));
        }
            
        
        
        //OK and now that the backbone atoms are in place, we can handle the sidechains and HA's
        for(int resNum=0; resNum<numRes; resNum++){
            //first translate into place...
            RigidBodyMotion motion = new RigidBodyMotion(new double[3], RotationMatrix.identity(), SCTranslations[resNum]);
            SidechainIdealizer.moveSidechain(residues.get(resNum), motion);
            //...now idealize
            SidechainIdealizer.idealizeSidechain(EnvironmentVars.resTemplates, residues.get(resNum));
            //and get gen chi1 to where it was before, so this BB motion commutes w/ sidechain dihedral changes
            GenChi1Calc.setGenChi1(residues.get(resNum), genChi1[resNum]);
        }
        
        curFreeDOFVals = x.toArray();
    }
    
    
    public ArrayList<BBFreeDOF> getDOFs(){
        return freeDOFs;
    }

    public double[][] getFreeDOFVoxel() {
        return freeDOFVoxel;
    }

    public List<Residue> getResidues() {
        return residues;
    }
    
    
    

    
    
}
