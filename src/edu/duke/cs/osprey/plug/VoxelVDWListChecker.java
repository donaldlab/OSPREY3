/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug;

import cern.colt.matrix.linalg.SingularValueDecomposition;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import edu.duke.cs.osprey.bbfree.BBFreeBlock;
import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import static edu.duke.cs.osprey.plug.LPChecks.canAddConstr;
import static edu.duke.cs.osprey.plug.VoxelVDWDistExplorer.getVDWRadius;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.Relationship;

/**
 *
 * For a specified voxel (bounds on specified dihedrals),
 * check if VDW interactions between the specified atoms
 * can be made to fall in range
 * FOR NOW use linear programming, (relaxing clashes??  or linearize if need, thus cutting vox)
 * 
 * @author mhall44
 */
public class VoxelVDWListChecker {
    
    protected static class DOFInterval {
        DegreeOfFreedom dof;
        double lb;
        double range;
        
        protected DOFInterval(DegreeOfFreedom dof, double lb, double range) {
            this.dof = dof;
            this.lb = lb;
            this.range = range;
        }
    }
    
    ArrayList<DOFInterval> dofIntervals = new ArrayList<>();
    ArrayList<Atom[]> interactingAtoms = new ArrayList<>();
    ArrayList<Residue> knownConfRes = null;//residues of known conf, for use in picking uc
    HashSet<Residue> shellResidues = null;
    
    ArrayList<LinearConstraint> voxLinConstr = new ArrayList<>();//non-box constr in RC definition
    
    
    static boolean onlyLC = true;
    //In 1CC8, Phe 9 HD2 has no contacts. 
    //This is a very well packed structure.  CD2 has contacts
    //So "realistic" contact enforcement should probs go to unified-atom level if that
    //w/ so many options it's hard to pick a uc specifically,
    //though could take the best and enforce that, or just enforce that there
    //be something available in the voxel
    //specific uc w/i vox less important than lc bc prob doesn't change much ruggedness-wise
    //contact lists for voxels admits more flexible rules (1 exception, etc.) while still efficient

    public VoxelVDWListChecker(ArrayList<DOFInterval> dofIntervals, ArrayList<Atom[]> interactingAtoms, 
            ArrayList<Residue> knownConfRes, ArrayList<Residue> shellRes, 
            ArrayList<LinearConstraint> voxLinConstr) {
        this.dofIntervals = dofIntervals;
        this.interactingAtoms = interactingAtoms;
        this.knownConfRes = knownConfRes;
        this.voxLinConstr = voxLinConstr;
        
        if(shellRes==null){
            int quack =5;
        }
        
        shellResidues = new HashSet(shellRes);
    }
    
    
    public VoxelVDWListChecker(){}
    
    //for piece-by-piece construction
    public final void addDOFInterval(DegreeOfFreedom dof, double lb, double range){
        dofIntervals.add(new DOFInterval(dof,lb,range));
    }
    
    public final void addAtomPair(Atom atom1, Atom atom2){
        interactingAtoms.add(new Atom[] {atom1,atom2});
    }
    
    public final void addVoxLinConstr(ArrayList<LinearConstraint> linConstr){
        voxLinConstr.addAll(linConstr);
    }

    //ALSO CONSTRUCTOR FROM RCs???  (but need to specify atoms somehow)
    //(will need calc of what atoms are w/i the VDW range+0.25 somewhere in the voxtopher)
    
    int numDOFs(){
        return dofIntervals.size();
    }
    
    boolean checkFeasibility(){
        return calcFeasiblePolytope(null)!=null;
    }
    
    
    double atomPairOverlap(Atom at1, Atom at2){
        double dist = VectorAlgebra.distance(at1.getCoords(), at2.getCoords());
        return getVDWRadius(at1) + getVDWRadius(at2) - dist;
    }
    
    
    ArrayList<LinearConstraint> calcFeasiblePolytope(ArrayList<String> atomPairNames){
        //If atomPairNames isn't null then record names of interacting atoms that generate inequalities there
        
        //first add lc
        ArrayList<LinearConstraint> curPolygon = voxelPolygon();

        if(atomPairNames!=null) {
            for (int n = 0; n < curPolygon.size(); n++)//non-steric constraints
                atomPairNames.add("VOXEL CONSTRAINT");
        }
                
        int numVDWPairs = interactingAtoms.size();
        for(int p=0; p<numVDWPairs; p++){            
            LinearConstraint lc = tooCloseConstr(interactingAtoms.get(p));
            
            
            boolean validLC = ! (lc.getCoefficients().getNorm()<1e-10);
            //DEBUG!!!!!  it's sil to penalize unchangeable bb clashes etc
            //and if they're always true constraints no point either
            //ACTUALLY NO THIS ALLOWS ALL KINDS OF RIDICULOUS ROTS
            //FOR NOW LET'S KILL CONF IF LC ALWAYS BAD (BELOW)
            //ONLY TOLERATING CLASHES INVOLVING ONLY UNMOVED ATOMS
            //(DOES NOT APPLY IN CATS RES)
            
            if(validLC){
                if(!canAddConstr(lc,curPolygon)){
                    return null;
                }

                curPolygon.add(lc);
                if(atomPairNames!=null)
                    atomPairNames.add(makeAtomPairNames(interactingAtoms.get(p)));
            }
            else if( (lc.getRelationship()==Relationship.GEQ) == (lc.getValue()>0) ){//constraint unsatisfiable
                if(!atomsFixed(interactingAtoms.get(p)))
                    return null;//KILL
            }
        }
        
        if(onlyLC)
            return curPolygon;
        
        //now pick the closest atoms to add uc
        //assuming right residue types already in place
        for(DOFInterval di : dofIntervals){
            di.dof.apply(di.lb+di.range/2);
        }
        
        HashMap<Atom,Atom> closestContactingAtom = new HashMap<>();
        HashMap<Atom,Double> biggestOverlap = new HashMap<>();
        HashMap<Atom,ArrayList<Atom>> partners = new HashMap<>();
        
        for(Atom[] atomPair : interactingAtoms){
            //possibly check if lc ok here before putting in
            
            if(!partners.containsKey(atomPair[0]))
                partners.put(atomPair[0], new ArrayList<>());
            partners.get(atomPair[0]).add(atomPair[1]);
            
            double overlap = atomPairOverlap(atomPair[0], atomPair[1]);
            if(biggestOverlap.containsKey(atomPair[0])){
                if(overlap<=biggestOverlap.get(atomPair[0]))
                    continue;//this overlap is not the biggest
            }
            biggestOverlap.put(atomPair[0], overlap);
            closestContactingAtom.put(atomPair[0], atomPair[1]);
        }
        
        //ok now for each atom in cur res that has contacts, see if it needs to have a good-vdw constraint
        for(Atom at1 : partners.keySet()){
            if(at1.bonds.size()>3)//Won't require VDW...4 bonds already "hem in" this atom
                continue;

            double closestOverlap = biggestOverlap.get(at1);
            boolean keepUC = true;
            if(closestOverlap<0){//there's some room...see if something else could get closer
                boolean hasCloseAtoms = false;
                double dist1 = getVDWRadius(at1) + 1 - closestOverlap;//1 is minimum VDW rad of an atom
                //if can get in another atom this close then we won't enforce the uc because it 
                //might not really be the closest
                for(Atom at2 : partners.get(at1)){//partners are close enough to block this other atom
                    double dist2 = getVDWRadius(at2) + 1;//can't put the other atom too close to a partners
                    ArrayList<double[]> borderSamps = sampleSphereIntersection(at1,dist1,at2,dist2);
                    if(!borderSamps.isEmpty())
                        hasCloseAtoms = true;
                    for(double[] samp : borderSamps){
                        if(!coordWithinSomeAtomVDWRad(samp,1)){//DEBUG!!! this is a very slow way to do this
                            //note we are assuming that coord will not overlap with any other atoms
                            //not just that nucleus is outside anything else's vdw rad
                            keepUC = false;
                            break;
                        }
                    }
                    if(!keepUC)
                        break;
                }
            
                if(!hasCloseAtoms)//no atoms close enough to at 1 to be a meaningful constr
                    keepUC = false;
            }
                
            if(keepUC){
                LinearConstraint uc = tooFarConstr(new Atom[] {at1,closestContactingAtom.get(at1)});
                boolean validUC = ! (uc.getCoefficients().getNorm()<1e-10);
                if(validUC){
                    if(!canAddConstr(uc,curPolygon)){
                        return null;
                    }
                    curPolygon.add(uc);//DEBUG!!!  since don't intend to use uc not including atom pair names here
                }
            }
        }
        
        //NOTE AT END OF ASSIGNMENTS, AND IDEALLY EVEN BEFORE THEN, WE CAN CHECK THAT EACH ATOM THAT SHOULD HAVE AN ATOM DOES HAVE ONE
        
        //if we got here we have a feasible pt and its containing polygon!!!
        return curPolygon;
    }
    
    
    String makeAtomPairNames(Atom[] atoms){
        //record names for atom pair
        return makeAtomName(atoms[0])+" ; "+makeAtomName(atoms[1]);
    }
    
    String makeAtomName(Atom atom){
        return atom.res.fullName+" , "+atom.name;
    }
    
    
    boolean atomsFixed(Atom[] pair){
        //DEBUG!!!!  Hacky way to deal with atoms that can't move out of a clash
        //so we don't penalize the clash
        for(DOFInterval di : dofIntervals){
            if(di.dof instanceof BBFreeDOF)//CATS can move everything in principle
                return false;
        }
        for(Atom at : pair){
            if(HardCodedResidueInfo.possibleBBAtomsLookup.contains(at.name))
                continue;
            if(at.name.equalsIgnoreCase("CB"))
                continue;//not technically backbone but can't move either
            if(shellResidues.contains(at.res))
                continue;
            return false;//atom can move 
        }
        return true;//didn't find a way for the atom to move
    }
    
    boolean coordWithinSomeAtomVDWRad(double coord[], double buffer){
        //DEBUG!!!!  SLO SLO SLO
        for(Residue res : knownConfRes){
            for(Atom at : res.atoms){
                double dist = VectorAlgebra.distance(at.getCoords(), coord);
                if(dist<getVDWRadius(at)+buffer)
                    return true;
            }
        }
        return false;
    }
    
    
    
    static ArrayList<double[]> sampleSphereIntersection(Atom at1,double rad1,Atom at2,double rad2){
        //sample point on the intersection of the spheres with specified radii
        //centered at the specified atoms
        ArrayList<double[]> ans = new ArrayList<>();
        double[] aaVec = VectorAlgebra.subtract(at2.getCoords(), at1.getCoords());
        double dist = VectorAlgebra.norm(aaVec);
        double cosAng1 = (rad1*rad1+dist*dist-rad2*rad2) / (2*rad1*dist);
        
        if(cosAng1>1 || cosAng1<-1)//spheres don't intersect
            return ans;
        
        double z = cosAng1*rad1;
        double r = rad1*Math.sqrt(1.-cosAng1*cosAng1);//radius of intersection circle
        double cen[] = VectorAlgebra.add(at1.getCoords(), VectorAlgebra.scale(aaVec,z/dist));//its center
        
        double x[] = VectorAlgebra.getPerpendicular(aaVec);
        VectorAlgebra.normalizeInPlace(x);
        double y[] = VectorAlgebra.cross(aaVec, x);
        VectorAlgebra.normalizeInPlace(y);
        
        double numSamps = 8;//DEBUG!!
        for(int s=0; s<numSamps; s++){
            double ang = s*2*Math.PI/numSamps;
            double[] samp = VectorAlgebra.scale(x, r*Math.cos(ang));
            VectorAlgebra.addInPlace(samp, VectorAlgebra.scale(y, r*Math.sin(ang)) );
            VectorAlgebra.addInPlace(samp, cen);
            ans.add(samp);
            
            //double checkr1 = VectorAlgebra.distance(samp, at1.getCoords());
            //double checkr2 = VectorAlgebra.distance(samp, at2.getCoords());
        }
        
        return ans;
    }
    
    
    //OK attempting to do UC by selecting with LC
    ArrayList<LinearConstraint> calcFeasiblePolytope3(){
        
        //first get lc polygon, to prune bad uc
        ArrayList<LinearConstraint> lcPolygon = voxelPolygon();
                
        int numVDWPairs = interactingAtoms.size();
        for(int p=0; p<numVDWPairs; p++){            
            LinearConstraint lc = tooCloseConstr(interactingAtoms.get(p));
            boolean validLC = ! (lc.getCoefficients().getNorm()<1e-10);
            
            if(validLC){
                if(!canAddConstr(lc,lcPolygon)){
                    return null;
                }
                lcPolygon.add(lc);
            }
        }
            
        
        //ok now get the real polygon
        ArrayList<LinearConstraint> curPolygon = (ArrayList<LinearConstraint>)lcPolygon.clone();
        for(int p=0; p<numVDWPairs; p++){            
            LinearConstraint uc = tooFarConstr(interactingAtoms.get(p));
            
            boolean validUC = ! (uc.getCoefficients().getNorm()<1e-10);
            
            if(validUC){
                if(canAddConstr(uc,lcPolygon)){
                    if(!canAddConstr(uc,curPolygon)){
                        return null;
                    }
                }
                curPolygon.add(uc);
            }
        }
        
        //if we got here we have a feasible pt and its containing polygon!!!
        return curPolygon;
    }
    
    
    ArrayList<LinearConstraint> calcFeasiblePolytope2(){
        //for geom search will want to return a polygon, and valid solns will be those
        //that are in the intersection of polygons.  
        
        //LET'S SEARCH THE VOXEL DEFINED AS THE WELL W/O CONVEX BLOCKAGES,
        //AS VIEWED FROM CENTER.  THUS LP SUFFICES, LINEARIZATION ALL FROM CENTER
        
        //ALTERNATIVE: DO THIS WITH CONVEX OPT
        //KEEP FEAS PT THE SAME IF POSSIBLE, ELSE MINIMIZE DIST TO VOX CENTER WI "POLYGON"
        //START W LINEAR CONSTR BUT SHOULDNT BE HARD TO ADD CONVEX QUADR IF NEEDED
        //ALBEIT NOT CLEAR IF WORTH IT SINCE CANT DO CONCAVE QUADR
        //AND THUS ARE BASICALLY REDEFINING VOXELS TO BE UNION OF CONVEX "SHADOW" OF FEASIBLE REGIONS
        //AS VIEWED FROM DIRECTION OF VOX CENTER
        
        //DoubleMatrix1D curFeasiblePt = voxelCenter();
        //DONT NEED TO MAINTAIN THIS IN LP/CENTER-BASED SHADOW FRAMEWORK
        
        
        ArrayList<LinearConstraint> curPolygon = voxelPolygon();
                
        int numVDWPairs = interactingAtoms.size();
        for(int p=0; p<numVDWPairs; p++){            
            LinearConstraint lc = tooCloseConstr(interactingAtoms.get(p));
            
            
            boolean validLC = ! (lc.getCoefficients().getNorm()<1e-10);
            //DEBUG!!!!!  it's sil to penalize unchangeable bb clashes etc
            //and if they're always true constraints no point either

            
            if(validLC){
                if(!canAddConstr(lc,curPolygon)){
                    return null;
                }
                /*if(!lc.isFeasible(curFeasiblePt))
                    curFeasiblePt = optimize(lc, curPolygon, curFeasiblePt);
                if(curFeasiblePt==null)
                    return false;*/

                curPolygon.add(lc);
            }
            //same for uc...
            //DEBUG!!!  See what happens if don't enforce uc
            /*LinearConstraint uc = tooFarConstr(interactingAtoms.get(p));
            if(!canAddConstr(uc,curPolygon))
                return null;
            
            
            curPolygon.add(uc);
                    */
            //OK empirically each res is always hemmed in by a wall of atoms
            //want contact with that wall but not clash (almost each possible contact with the wall should occur--empirically they do)
            //For local pruning it's essential to choose wall atoms by local reasoning!  (w/i pair/tuple)
            //Can enforce lc for any atoms that may come close
            //For uc, could try enforcing any constraints compatible with lc that are possible
            //(maybe with buffer) w/i voxel
            //But best way is likely to only enforce constr between atom and closest atom in other res
            //or rest of res tuple (can put dummy atoms around boundary of interface,
            //at closest distance possible, to rule out atom not being closest)
            //dummy atoms go away as reach full conf
            //could even try Delaunay, enforce neighbors
            //this analysis can be done at central conf (really shouldn't qualitatively change much across conf)
            //or at a couple samples
            //EG in Met 13 theres trouble just from chi1 bc wants to enforce constr w/ multiple atoms
            //in Leu 40 cd methyl, and these end up being contradictory, but 
            //going by distance (ADJUSTED FOR VDWRAD) the cd vdw clearly dominates
        }
        
        //if we got here we have a feasible pt and its containing polygon!!!
        return curPolygon;
    }
        
    

    
    double[] DOFsAtAtomPairDist(Atom[] atomPair, double targetDist, double[] initGuess){
        //Starting at initGuess, find values for the DOFs that put the specified atoms
        //the given target distance apart
        double[] x = initGuess.clone();
        double initDist = calcAtomPairDist(x,atomPair);
        double curDist = initDist;
        
        double w[] = new double[numDOFs()];
        for(int dof=0; dof<numDOFs(); dof++)
            w[dof] = dofIntervals.get(dof).range;
        
        while( Math.abs(targetDist-curDist)  > 0.01 ){
            double grad[] = calcAtomPairDistGrad(x,atomPair);
            
            double derivsq = 0;
            for(int dof=0; dof<numDOFs(); dof++)
                //derivsq += grad[dof]*grad[dof];
                derivsq += w[dof]*w[dof]*grad[dof]*grad[dof];
            
            if(derivsq==0)//DOFs don't move dist, so can't reach target by moving DOFs
                return null;
            
            for(int dof=0; dof<numDOFs(); dof++){
                //x[dof] -= (curDist-targetDist) * grad[dof] / derivsq;
                x[dof] -= (curDist-targetDist) * grad[dof]*w[dof]*w[dof] / derivsq;
                
                DOFInterval di = dofIntervals.get(dof);
                if(x[dof]<di.lb || x[dof]>di.lb+di.range){
                    //DEBUG!!!  if we get outside the voxel probably we can't reach this boundary
                    return null;
                }
            }
            
            curDist = calcAtomPairDist(x,atomPair);
        }
        return x;
    }
    
    double[] DOFsAtAtomPairDist(Atom[] atomPair, double targetDist){
        //Start at voxel center by default
        return DOFsAtAtomPairDist(atomPair, targetDist, voxelCenter());
    }
    
    LinearConstraint linearizedAtomDistConstr(Atom[] atomPair, double targetDist, boolean lbConstr){
        //Linearized constraint on atom-atom distance
        
        double cen[] = voxelCenter();
        double x[] = DOFsAtAtomPairDist(atomPair, targetDist, cen);
        
        if(x==null){
            //couldn't get to boundary, so let's return an always-true or always-false
            double initDist = calcAtomPairDist(cen, atomPair);
            if(lbConstr==(initDist>targetDist))//constraint satisfied in voxel (seemingly by whole vox)
                return new LinearConstraint(new double[numDOFs()], Relationship.GEQ, -1.);//always true
            else
                return new LinearConstraint(new double[numDOFs()], Relationship.GEQ, 1.);//always false
        }
        
        double curDist = calcAtomPairDist(x, atomPair);
        double grad[] = calcAtomPairDistGrad(x,atomPair);//NUM FOR NOW
        double targ = targetDist-curDist;
        for(int dof=0; dof<numDOFs(); dof++)
            targ += grad[dof]*x[dof];
        
        //OK we should now be close enough to get a good approximate constr
        if(lbConstr)
            return new LinearConstraint(grad, Relationship.GEQ, targ);
        else//we are upper-bounding the distance
            return new LinearConstraint(grad, Relationship.LEQ, targ);
    }
    
    
    double calcAtomPairDist(double[] DOFVals, Atom[] atomPair){
        for(int dofNum=0; dofNum<numDOFs(); dofNum++)
            dofIntervals.get(dofNum).dof.apply( DOFVals[dofNum] );
        
        return VectorAlgebra.distance(atomPair[0].getCoords(), atomPair[1].getCoords());
    }
    
    double[] calcAtomPairDistGrad(double[] DOFVals, Atom[] atomPair){
        double step = 1e-6;
        int numDOFs = numDOFs();
        double ans[] = new double[numDOFs];
        for(int dof=0; dof<numDOFs; dof++){
            DOFVals[dof] += step;
            double upVal = calcAtomPairDist(DOFVals,atomPair);
            DOFVals[dof] -= 2*step;
            double downVal = calcAtomPairDist(DOFVals,atomPair);
            DOFVals[dof] += step;
            ans[dof] = (upVal-downVal)/(2*step);
        }
        return ans;
    }
    
    
    
    LinearConstraint tooCloseConstr(Atom[] atomPair){
        double idealDist = getVDWRadius(atomPair[0]) + getVDWRadius(atomPair[1]);
        double buffer = getClashBuffer(atomPair[0], atomPair[1]);
        return linearizedAtomDistConstr(atomPair, idealDist-buffer, true);
    }
    
    public static double getTargetDist(Atom atom1, Atom atom2){
        //how close can atom1, atom2 be without clashing
        return getVDWRadius(atom1)+getVDWRadius(atom2)-getClashBuffer(atom1,atom2);
    }
    
    static double getClashBuffer(Atom atom1, Atom atom2){
        //figure out the amount of VDW overlap allowed for the two atoms without 
        //begin considered a clash
        //follows Richardsons' "Contact-dot Surfaces with Explicit H Atoms" paper
        //normally 0.4 A, but lower for hydrogen bonds (specified as O, His N, S, or aromatic ring
        //face for acceptor, CH excluded as donor so must be NH or OH)
        //overlap is 0.6 A for most hydrogen bonds, 0.8 A for salt bridge
        if(atom1.elementNumber==1){//switch for easier analysis
            Atom tmp = atom1;
            atom1=atom2;
            atom2=tmp;
        }
        
        if(atom2.elementNumber==1){
            int HBondedElem = atom2.bonds.get(0).elementNumber;
            if(HBondedElem==7 || HBondedElem==8){
                if(atom1.elementNumber==8||atom1.elementNumber==16){                    
                    //OK we now have a hydrogen bond
                    //DEBUG!!! should probably count aromatic C or HIS N too
                    //and cover all salt bridges--here just ammonium or arginine to carboxylate
                    //full version should probably pre-go through molecule
                    //and assign each atom a probe radius, H-bond eligibility, and charged status
                    if(isCarboxylateO(atom1) && (isAmmoniumH(atom2)||isArgChargedH(atom2)))
                        return 0.8;//salt bridge
                    else
                        return 0.6;
                }
            }
        }
        
        return 0.4;
    }
    
    
    private static boolean isCarboxylateO(Atom atom){
        if(atom.elementNumber!=8)
            return false;
        if(atom.bonds.size()!=1)
            return false;
        Atom C = atom.bonds.get(0);
        if(C.elementNumber!=6)
            return false;
        if(C.bonds.size()!=3)
            return false;
        int countO = 0;
        for(Atom atom2 : C.bonds){
            if(atom2.bonds.size()==1 && atom2.elementNumber==8)
                countO++;
        }
        return countO==2;
    }
    
    
    private static boolean isAmmoniumH(Atom atom){
        if(atom.elementNumber!=1)
            return false;
        Atom N = atom.bonds.get(0);
        if(N.elementNumber!=7)
            return false;
        return N.bonds.size()==4;
    }
    
    private static boolean isArgChargedH(Atom atom){
        if(!atom.res.template.name.equalsIgnoreCase("ARG"))
            return false;
        if(atom.elementNumber!=1)
            return false;
        return atom.name.equalsIgnoreCase("HE") || atom.name.startsWith("HH");
    }
    
    
    LinearConstraint tooFarConstr(Atom[] atomPair){
        double idealDist = getVDWRadius(atomPair[0]) + getVDWRadius(atomPair[1]);
        return linearizedAtomDistConstr(atomPair, idealDist+0.25, false);
    }
    

   
    
    ArrayList<LinearConstraint> voxelPolygon(){
        ArrayList<LinearConstraint> ans = new ArrayList<>();
        for(int dof=0; dof<numDOFs(); dof++){
            double unitVec[] = new double[numDOFs()];
            unitVec[dof] = 1;
            DOFInterval di = dofIntervals.get(dof);
            ans.add(new LinearConstraint(unitVec,Relationship.GEQ,di.lb));
            ans.add(new LinearConstraint(unitVec,Relationship.LEQ,di.lb+di.range));
        }
        ans.addAll(voxLinConstr);
        return ans;
    }
    
    double[] voxelCenter(){
        double ans[] = new double[numDOFs()];
        HashMap<DegreeOfFreedom,Double> specialCenter = calcSpecialCenter();
        for(int dof=0; dof<numDOFs(); dof++){
            DOFInterval di = dofIntervals.get(dof);
            if(specialCenter.containsKey(di.dof))
                ans[dof] = specialCenter.get(di.dof);
            else
                ans[dof] = di.lb + di.range/2;
        }
        
        return ans;
    }
    
    
    public HashMap<DegreeOfFreedom,Double> calcSpecialCenter(){
        //center for some degrees of freedom may, due to non-box constr, 
        //differ from center of box constr
        //NOTE: the map only contains those degrees of freedom where center is special
        //(not just middle of lb,ub)
        //DEBUG!!! This assumes the particular form of lin constr (pairs of parallel constr,
        //DOFs not in constrained space assumed to center at 0) used in basis big CATS stuff
        
        HashMap<DegreeOfFreedom,Double> ans = new HashMap<>();
        int numSpecialDOFConstr = voxLinConstr.size()/2;
        if(numSpecialDOFConstr==0)
            return ans;
        
        HashSet<Integer> specialDOFSet = new HashSet<>();
        for(int spesh=0; spesh<numSpecialDOFConstr; spesh++){
            RealVector coeff = voxLinConstr.get(2*spesh).getCoefficients();
            if(coeff.getDistance(voxLinConstr.get(2*spesh+1).getCoefficients()) > 0)
                throw new RuntimeException("ERROR: voxLinConstr not paired like expected");
            
            for(int d=0; d<dofIntervals.size(); d++){
                if(coeff.getEntry(d)!=0)
                    specialDOFSet.add(d);
            }
        }
        
        int numSpecialDOFs = specialDOFSet.size();
        ArrayList<Integer> specialDOFList = new ArrayList<>();
        specialDOFList.addAll(specialDOFSet);
        
        DoubleMatrix2D specialConstr = DoubleFactory2D.dense.make(numSpecialDOFConstr, numSpecialDOFs);
        //we will constrain values first of the lin combs of DOFs given in voxLinConstr,
        //then of lin combs orth to those
        DoubleMatrix2D target = DoubleFactory2D.dense.make(numSpecialDOFs,1);
        for(int spesh=0; spesh<numSpecialDOFConstr; spesh++){
            RealVector coeffs = voxLinConstr.get(2*spesh).getCoefficients();
            for(int d=0; d<numSpecialDOFs; d++)
                specialConstr.set(spesh, d, coeffs.getEntry(specialDOFList.get(d)));
            target.set( spesh, 0, 0.5*(voxLinConstr.get(2*spesh).getValue()+voxLinConstr.get(2*spesh+1).getValue()) );
        }
        
        //remaining targets (orth components to voxLinConstr) are 0
        DoubleMatrix2D orthComponents = getOrthogVectors(specialConstr);
        DoubleMatrix2D M = DoubleFactory2D.dense.compose(
                new DoubleMatrix2D[][] {
                    new DoubleMatrix2D[] {specialConstr},
                    new DoubleMatrix2D[] {Algebra.DEFAULT.transpose(orthComponents)}
                }
            );
        
        DoubleMatrix1D specialDOFVals = Algebra.DEFAULT.solve(M, target).viewColumn(0);
        
        for(int d=0; d<numSpecialDOFs; d++){
            ans.put(dofIntervals.get(specialDOFList.get(d)).dof, specialDOFVals.get(d));
        }
        
        return ans;
    }


    public static DoubleMatrix2D getOrthogVectors(DoubleMatrix2D M){
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
    
}
