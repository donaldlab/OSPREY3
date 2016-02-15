/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import edu.duke.cs.osprey.bbfree.BBFreeBlock;
import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.EllipseCoordDOF;
import edu.duke.cs.osprey.dof.MoveableStrand;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.dof.deeper.perts.Perturbation;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.*;
import java.io.Serializable;
import java.util.ArrayList;


/**
 *
 * @author mhall44
 */
public class PositionConfSpace implements Serializable {
    //This class defines the conformational space of a flexible residue
    //including allowed amino-acid types, and rotamers/RCs for each type 
    //subclass PositionConfSpace to make super-residues with super-RCs...
    
    
    public ArrayList<RC> RCs = new ArrayList<>();
    
    public Residue res;//The residue involved
    public int designIndex;
    
    static double dihedFlexInterval = 9;// +/- 9 degree sidechain dihedral continuous flexibility...
    //later can allow this to vary across different dihedrals
    
    static double ellipseAngMax = Math.PI;
    static double ellipseFinAngMax = Math.PI * 2;
    static double ellipseMin = 0;
    
    ArrayList<EllipseCoordDOF> ellipsoidalDOFs = new ArrayList<>();
    
    public PositionConfSpace(int pos, Residue res, ArrayList<DegreeOfFreedom> resDOFs, ArrayList<String> allowedAAs, 
            boolean contSCFlex, ArrayList<DegreeOfFreedom> strandDOFs, 
            ArrayList<Perturbation> perts, ArrayList<ArrayList<double[]>> pertIntervals, 
            ArrayList<ArrayList<int[]>> pertStates, BBFreeBlock bfb, boolean useEllipses ) {
        
        //We'll start with just one RC for each rotamer
        //But in general there are a lot of options for RCs...
        
        ResidueTemplateLibrary templateLib = EnvironmentVars.resTemplates;
        this.res = res;
        designIndex = pos;
        
        
        if(pertStates==null){//no DEEPer flexibility...
            pertStates = new ArrayList<>();
            pertStates.add(null);
        }
        
        
        for(String AAType : allowedAAs){
            int numDihedrals = templateLib.numDihedralsForResType(AAType);
            
        	// PGC 2015: Support backbone dependent rotamers.  
        	//	Compute phi and psi, necessary for backbone dependent rotamers.        
        	double phipsi [] = Protractor.getPhiPsi(this.res);
        	int numRot = templateLib.numRotForResType(designIndex, AAType, phipsi[0], phipsi[1]);
            
            //resDOFs is all sidechain DOFs, for now
            ArrayList<DegreeOfFreedom> dofListForRot = new ArrayList<>();
            for(int dih=0; dih<numDihedrals; dih++)//get the first numDihedrals dihedrals
                dofListForRot.add(resDOFs.get(dih));
                       
            
            for(ArrayList<int[]> pertState : pertStates){
                
                if(AAType.equalsIgnoreCase("PRO")){//special case: has ring pucker
                    //If PRO is set to be flexible we'll assume this includes pucker flexibility
                    //(the main flexibility of proline)
                    for( int puckerVal : new int[]{0,1} ){//create both puckers
                        createRC(0, AAType, -1, contSCFlex, dofListForRot, puckerVal,
                                strandDOFs, bfb, pertState, perts, pertIntervals, false);
                    }
                }
                
                else if(numRot==0){//ALA or GLY: no rotamers or dihedrals, so create a single rigid RC (or one for each pert state)
                    createRC(0, AAType, -1, contSCFlex, dofListForRot, -1, strandDOFs, bfb, 
                            pertState, perts, pertIntervals, false);
                }
                
                else {
                    boolean AATypeUsesEllipses = useEllipses && (numDihedrals>1);//ellipses only meaningful if > 1 dihedral
                    
                    for(int rot=0; rot<numRot; rot++){
                        createRC(numDihedrals, AAType, rot, contSCFlex, dofListForRot, -1, strandDOFs, 
                                bfb, pertState, perts, pertIntervals, AATypeUsesEllipses);
                    }
                }
            }
        }
        
    }
    
    
    private void createRC(int numDihedrals, String AAType, int rot, boolean contSCFlex, ArrayList<DegreeOfFreedom> dofListForRot,
            int proPucker, ArrayList<DegreeOfFreedom> strandDOFs, BBFreeBlock bfb, ArrayList<int[]> pertState, 
            ArrayList<Perturbation> perts, ArrayList<ArrayList<double[]>> pertIntervals, boolean useEllipses){
        
        //Create an RC with the specified rotamer and (if set) perturbation state
        
        ResidueTemplateLibrary templateLib = EnvironmentVars.resTemplates;
        //create RC
        ArrayList<Double> dofLB = new ArrayList<>();//lower bounds on each DOF for this RC
        ArrayList<Double> dofUB = new ArrayList<>();//upper bounds
        
        
    	// PGC 2015: Support backbone dependent rotamers.  
    	//	Compute phi and psi, necessary for backbone dependent rotamers.        
    	double phipsi [] = Protractor.getPhiPsi(this.res);
        double dihValues[] = new double[numDihedrals];
        //we'll start with the sidechain dihedral DOFs

        for(int dih=0; dih<numDihedrals; dih++){
            dihValues[dih] = templateLib.getDihedralForRotamer(designIndex, AAType,phipsi[0], phipsi[1], rot,dih);
        }
        
        ArrayList<EllipseCoordDOF> ellCoords = makeEllCoords(useEllipses, dihValues, dofListForRot);
        
        boundRotDOFs(dofLB, dofUB, numDihedrals, contSCFlex, dihValues, ellCoords, useEllipses);
        
        ArrayList<DegreeOfFreedom> dofListForRC = new ArrayList<>();
        if(useEllipses){
            dofListForRC.addAll(ellCoords);
            ellipsoidalDOFs.addAll(ellCoords);
        }
        else {   
            dofListForRC.addAll(dofListForRot);
        }
        
        //Put in proline pucker if this is a proline
        if(AAType.equalsIgnoreCase("PRO")){
            dofListForRC.add(res.pucker);
            dofLB.add( (double) proPucker );
            dofUB.add( (double) proPucker );
        }
        
        addStrandDOFs(strandDOFs, dofListForRC, dofLB, dofUB);
        addBFBDOFs(bfb, dofListForRC, dofLB, dofUB);
        addDEEPerDOFs(pertState, perts, pertIntervals, dofListForRC, dofLB, dofUB);
        
        RC newRC = new RC(AAType, rot, dofListForRC, dofLB, dofUB, RCs.size());
        
        RCs.add(newRC);
    }
    
    
    //methods to add DOFs besides sidechain flexibility to an RC
    
    private void addStrandDOFs(ArrayList<DegreeOfFreedom> strandDOFs, ArrayList<DegreeOfFreedom> dofListForRC,
            ArrayList<Double> dofLB, ArrayList<Double> dofUB ){
        //Add strand rotation/translation DOFs to an RC, if there are any, putting them in dofListForRC, dofLB, and dofUB
        //these are DOFs whose strand includes this residue!
        for(DegreeOfFreedom strandDOF : strandDOFs){
            dofListForRC.add(strandDOF);
            double[] bounds = MoveableStrand.getStrandDOFBounds(strandDOF);
            dofLB.add(bounds[0]);
            dofUB.add(bounds[1]);
        }
    }
    
    
    private void addBFBDOFs(BBFreeBlock bfb, ArrayList<DegreeOfFreedom> dofListForRC,
            ArrayList<Double> dofLB, ArrayList<Double> dofUB){
        //Add free-BB block DOFs to an RC, if any
        
        if(bfb!=null){//This res is in a free-backbone block
            ArrayList<BBFreeDOF> bbFreeDOFs = bfb.getDOFs();
            double freeDOFVoxel[][] = bfb.getFreeDOFVoxel();
            
            for(int dnum=0; dnum<bbFreeDOFs.size(); dnum++){
                dofListForRC.add(bbFreeDOFs.get(dnum));
                dofLB.add(freeDOFVoxel[0][dnum]);
                dofUB.add(freeDOFVoxel[1][dnum]);
            }
        }
    }
    
    
    private void addDEEPerDOFs(ArrayList<int[]> pertState, 
            ArrayList<Perturbation> perts, ArrayList<ArrayList<double[]>> pertIntervals, 
            ArrayList<DegreeOfFreedom> dofListForRC,
            ArrayList<Double> dofLB, ArrayList<Double> dofUB){
        //Add DEEPer DOFs to an RC, if any
        
        if(pertState != null) {
            //need to add DEEPer DOFs
            
            for(int[] singlePertState : pertState){
                int pertNum = singlePertState[0];//number (in perts) of the perturbation we're adding
                int pertStateNum = singlePertState[1];//state of this perturbation
                dofListForRC.add(perts.get(pertNum));
                
                double[] pertInterval = pertIntervals.get(pertNum).get(pertStateNum);//interval for this perturbation in this state
                dofLB.add(pertInterval[0]);
                dofUB.add(pertInterval[1]);
            }
        }
    }
    
    
    private void boundRotDOFs(ArrayList<Double> dofLB, ArrayList<Double> dofUB, int numDihedrals,
            boolean contSCFlex, double[] dihValues, ArrayList<EllipseCoordDOF> ellCoords, boolean useEllipses){
        //Put the bounds on a rotamer's degrees of freedom (dihedrals or ellipsoidal DOFs) in dofLB and dofUB
        
        for(int dih=0; dih<numDihedrals; dih++){
            if(useEllipses){
                if (contSCFlex) {
                    dofLB.add(ellipseMin);
                    double radMax = 30;
                    dofUB.add((dih==0) ? radMax : // TODO: make this better
                                    (dih==numDihedrals-1) ? ellipseFinAngMax : ellipseAngMax);
                } else {
                    dofLB.add(ellCoords.get(dih).getCurVal());
                    dofUB.add(ellCoords.get(dih).getCurVal());
                }
            }
            else {
                if(contSCFlex){//allow continuous flexibility up to dihedFlexInterval in each direction
                    dofLB.add(dihValues[dih]-dihedFlexInterval);
                    dofUB.add(dihValues[dih]+dihedFlexInterval);
                }
                else {
                    dofLB.add(dihValues[dih]);
                    dofUB.add(dihValues[dih]);
                }
            }
        }
    }
    
    
    
    private ArrayList<EllipseCoordDOF> makeEllCoords(boolean useEllipses, double dihValues[], 
            ArrayList<DegreeOfFreedom> dofListForRot){
        //build the ellipsoidal coordinates for an RC.  Empty list if not using ellipses.  
        
        ArrayList<EllipseCoordDOF> ellCoords = new ArrayList<>();

        if(useEllipses){
            // generate the list of ellipsoidal DOFs
            // TODO: move getellipsoidalcoords to ellipsetransform
            double[] ellValues = getEllipsoidalCoords(dihValues);
            DoubleMatrix2D A = DoubleFactory2D.dense.identity(ellValues.length);
            DoubleMatrix1D c = DoubleFactory1D.dense.make(new double[ellValues.length]);
            for (int i=0; i<ellValues.length; i++) {
                    EllipseCoordDOF ellDOF = new EllipseCoordDOF(
                                    (i==0),
                                    i,
                                    ellValues[i],
                                    A,
                                    dofListForRot,
                                    dihValues);
                    ellCoords.add(ellDOF);

            }
        }
        
        return ellCoords;
    }
    
    
    public double[] getEllipsoidalCoords(double[] dihedrals) {
    	
        if (dihedrals.length==0) { return new double[0]; }
    	
    	// for now we're just using the unit sphere
    	DoubleMatrix2D A = DoubleFactory2D.dense.identity(dihedrals.length);
    	DoubleMatrix1D c = DoubleFactory1D.dense.make(new double[dihedrals.length]);
    	
        EigenvalueDecomposition evd = new EigenvalueDecomposition(A);
        DoubleMatrix2D Q = evd.getV();
        DoubleMatrix2D L = evd.getD();
        DoubleMatrix2D qT = Q.viewDice().copy();
        Algebra alg = new Algebra();
    	
    	// first transform the cartesian coordinates based on the ellipse
    	double[] s = new double[dihedrals.length];
    	for (int i=0; i<dihedrals.length; i++) { 
    		s[i] = dihedrals[i]-c.get(i);
    	}
    	double[] u = alg.mult(qT, DoubleFactory1D.dense.make(s)).toArray();
    	double[] x = new double[u.length];
    	for (int i=0; i<u.length; i++) {
    		x[i] = u[i]/Math.sqrt(L.get(i, i));
    	}    	
    	dihedrals = x;
    	
    	// now get elliptical coordinates
/*    	double radius = 0;
    	for (double d : dihedrals) { radius += d*d; }
    	radius = Math.sqrt(radius);*/
    	int n = dihedrals.length;
    	double[] phi = new double[n-1];
    	for (int i=0; i<n-1; i++) {
    		double d=0;
    		for (int j=i; j<n; j++) { d += dihedrals[j]*dihedrals[j]; }
    		double quot = dihedrals[i]/Math.sqrt(d);
    		phi[i] = Math.acos(quot);
    	}
    	if (dihedrals[n-1] < 0) { phi[n-2] = 2*Math.PI - phi[n-2]; }
    	double[] ellCoords = new double[n];
    	ellCoords[0] = 0; //radius;
    	for (int i=1; i<n; i++) { ellCoords[i] = phi[i-1]; }
    	return ellCoords;
    }

    public ArrayList<EllipseCoordDOF> getEllipsoidalArray() {
    	return this.ellipsoidalDOFs;
    }

    
}
