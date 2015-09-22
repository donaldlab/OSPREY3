/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.bbfree.BBFreeBlock;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.EllipseCoordDOF;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.dof.MoveableStrand;
import edu.duke.cs.osprey.dof.ProlinePucker;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.dof.StrandRotation;
import edu.duke.cs.osprey.dof.StrandTranslation;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.dof.deeper.SidechainIdealizer;
import edu.duke.cs.osprey.dof.deeper.perts.Perturbation;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MolecEObjFunction;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.StringParsing;
import java.io.Serializable;
import java.util.ArrayList;
/**
 *
 * @author hmn5
 * @author mhall44    
 * This class represents the conformational search space for a design
    used for GMEC-based design, K*, or anything else we want to do with a conformational space
    This class just defines the conformational space itself (the molecule + all kinds of flexibility
    nd possible mutations, and how these are represented as RCs, etc.)
    This class can be put in an AnnotatedConfSpace to add annotations like what RCs are pruned,
    what their pairwise energies are, etc.  
 */
public class ConfSpaceSuper {//extends ConfSpace{
//This class implementation is going to copy, rather than truly extend ConfSpace because we will probably
//not need both. If we do want to keep both, just remove the redudancies between the two
 
    

    public Molecule m;
    //The molecule will be composed of residues. 
    //It will have one set of coordinates, which are stored in the residues to make mutation
    //and pairwise energy computation easy (no need for complicated partial arrays, subtracting off
    //template energies, etc., and mutation will be very quick and will only require
    //AMBER reinitialization for the affected residue)
    
    //If we need to keep a copy of, say, the original coordinates, we can have a second molecule origMolec
    //for that
    //Once loaded, the molecule can only be changed by functions overriding DegreeOfFreedom.applyValue
    //this will give us a uniform framework for applying conformational and sequence changes
    //and each DOF will store its value, so we can easily and unambiguously look up the state of the
    //molecule at any given time
   
    //list of conformational degrees of freedom
    public ArrayList<DegreeOfFreedom> confDOFs = new ArrayList<>();
  
    //list of sequence degrees of freedom (residue types for the mutable residues)
    public ArrayList<ResidueTypeDOF> mutDOFs = new ArrayList<>();
    
    
    
    public ArrayList<PositionConfSpaceSuper> posFlex = new ArrayList<>();
    //defines the flexible positions and what RCs they have
    //generally each position is a residue, but it could be more than one ("super-residue" with "super-RCs")

    
    public int numPos;//number of flexible positions
    
    public boolean useEllipses = false;    

    //Default constructor will simply create posFlex with one residue per position
    //The only difference is everything is embedded in a list now so that more than
    //One residue at each position could be added.
    public ConfSpaceSuper(
            String PDBFile, 
            ArrayList<String> flexibleRes, 
            ArrayList<ArrayList<String>> allowedAAs, 
            boolean addWT, 
            boolean contSCFlex,
            DEEPerSettings dset,
            ArrayList<String[]> moveableStrands,
            ArrayList<String[]> freeBBZones,
            boolean ellipses){

        this.useEllipses =  ellipses;
        this.numPos = flexibleRes.size();
        
        //read the structure and assign templates, deleting unassignable res...
        m = PDBFileReader.readPDBFile(PDBFile);

        //Make all the degrees of freedom
        //start with proline puckers (added to res)
        makeRingPuckers(allowedAAs, flexibleRes);
        
        ArrayList<ArrayList<DegreeOfFreedom>> singleResDOFs = new ArrayList<>();
        for(int pos = 0; pos < numPos; pos++){
            
            Residue res = m.getResByPDBResNumber(flexibleRes.get(pos));
            if(addWT){//at this point, m has all wild-type residues, so just see what res is now
                String wtName = res.template.name;
                if(!allowedAAs.get(pos).contains(wtName))
                    allowedAAs.get(pos).add(wtName);
            }


            ArrayList<DegreeOfFreedom> resDOFs = mutableResDOFs(res, allowedAAs.get(pos));//add mutation and dihedral confDOFs for residue
            
            ResidueTypeDOF resMutDOF = (ResidueTypeDOF)resDOFs.remove(0);//first mutable pos DOF is the mutation-type DOF
            mutDOFs.add(resMutDOF);
            
            singleResDOFs.add(resDOFs);
            
        }
        
        //now rigid-body strand motions...
        ArrayList<DegreeOfFreedom> strandDOFs = strandMotionDOFs(moveableStrands, flexibleRes);
        confDOFs.addAll(strandDOFs);
        
        //...and perturbations
        //standardize conformations first since we'll record initial resBBState here
        standardizeMutatableRes(allowedAAs, flexibleRes);

        ArrayList<Perturbation> perts = dset.makePerturbations(m);//will make pert block here
        confDOFs.addAll(perts);
        
        //DEBUG!!!!!!!
        //TRYING BFB ON ALL FLEX RES!!!
        /*ArrayList<Residue> bfbRes = new ArrayList<>();
        for(String fr : flexibleRes)
            bfbRes.add( m.getResByPDBResNumber(fr) );
        BBFreeBlock bfb = new BBFreeBlock(bfbRes);
        confDOFs.addAll( bfb.getDOFs() );*/
        //DEBUG!!!
        ArrayList<BBFreeBlock> bfbList = getBBFreeBlocks(freeBBZones,flexibleRes);
        for(BBFreeBlock bfb : bfbList){
            confDOFs.addAll( bfb.getDOFs() );        
        }
        
        //OK now make RCs using these DOFs
        for(int pos=0; pos<numPos; pos++){
            ArrayList<Residue> resList = new ArrayList<>();
            Residue res = m.getResByPDBResNumber( flexibleRes.get(pos) );
            //just add one residue per position initially
            resList.add(res);
            
            ArrayList<DegreeOfFreedom> resDOFs = singleResDOFs.get(pos);
            ArrayList<DegreeOfFreedom> resStrandDOFs = strandDOFsForRes(res,strandDOFs);//and check that all res flex on moving strands!
            
            BBFreeBlock curBFB = getCurBFB(bfbList,res);

            //HMN: Create allowedAA at position
            ArrayList<ArrayList<String>> allowedAAsAtPosition = new ArrayList<ArrayList<String>>();
            allowedAAsAtPosition.add(allowedAAs.get(pos)); 
            //HMN: Create PosDOFS (ResDOFs per position)
            ArrayList<ArrayList<DegreeOfFreedom>> posDOFs = new  ArrayList<>();
            posDOFs.add(resDOFs);
            //HMN: getPertStatesPos = getPerStates(resNum) for each resNum in position
            ArrayList<ArrayList<ArrayList<int[]>>> pertStatesPos = new ArrayList<>();
            pertStatesPos.add(dset.getPertStates(pos));
            //HMN: BFB for each res in pos
            ArrayList<BBFreeBlock> curBFBPerRes = new ArrayList<>();
            curBFBPerRes.add(getCurBFB(bfbList, res));
            
            ArrayList<Integer> DOFIndices = new ArrayList<>();
            DOFIndices.add(pos);
            
            PositionConfSpaceSuper rcs = new PositionConfSpaceSuper(resList, posDOFs, allowedAAsAtPosition, DOFIndices, contSCFlex,
                       resStrandDOFs, perts, dset.getPertIntervals(), pertStatesPos, curBFBPerRes, useEllipses);
            posFlex.add(rcs);
                        
            if (useEllipses) {
            	confDOFs.addAll(rcs.getEllipsoidalArray());
            } else {
            	confDOFs.addAll(resDOFs);
            }
            
        }        
    }
        private ArrayList<BBFreeBlock> getBBFreeBlocks(ArrayList<String[]> freeBBZones, ArrayList<String> flexibleRes){
        //create a BFB for each (start res, end res) pair.  PDB residue numbers provided.  
        ArrayList<BBFreeBlock> ans = new ArrayList<>();
        
        for(String[] termini : freeBBZones){
            ArrayList<Residue> curBFBRes = resListFromTermini(termini, flexibleRes);
            BBFreeBlock bfb = new BBFreeBlock(curBFBRes);
            ans.add(bfb);
        }
        
        return ans;
    }
    
    private BBFreeBlock getCurBFB(ArrayList<BBFreeBlock> bfbList, Residue res){
        //If res is in one of the BFB's, return that BFB; else return null;
        for(BBFreeBlock bfb : bfbList){
            if(bfb.getResidues().contains(res)){//current res is in this BFB
                return bfb;
            }
        }
        
        return null;//not in a BFB
    }
    
    
    private void makeRingPuckers(ArrayList<ArrayList<String>> allowedAAs, ArrayList<String> flexibleRes){
        //For each residue, create a Pro pucker DOF if Pro is allowed
        
        for(int resNum=0; resNum<numPos; resNum++){
            Residue res = m.getResByPDBResNumber(flexibleRes.get(resNum) );
            
            for(String AAType : allowedAAs.get(resNum)){
                if(AAType.equalsIgnoreCase("PRO")){
                    res.pucker = new ProlinePucker(res);
                    break;
                }
            }
        }
    }
    
    private ArrayList<Residue> resListFromTermini(String[] termini, ArrayList<String> flexibleRes) {
        //Return a list of residues given the PDB numbers of the first and last
        //All of these residues are expected to be flexible (used for rot/trans strands and BBFreeBlocks)

        ArrayList<Residue> resList = new ArrayList<>();//res in current moving strand

        Residue curRes = m.getResByPDBResNumber(termini[0]);
        resList.add(curRes);

        while (!curRes.getPDBResNumber().equalsIgnoreCase(termini[1])) {//not at other end

            int curIndex = curRes.indexInMolecule;
            if (curIndex == m.residues.size() - 1) {
                throw new RuntimeException("ERROR: Reached end of molecule"
                        + " in rot/trans strand or BBFreeBlock without finding res " + termini[1]);
            }

            curRes = m.residues.get(curRes.indexInMolecule + 1);
            String curPDBNum = curRes.getPDBResNumber();
            if (!flexibleRes.contains(curPDBNum)) {
                throw new RuntimeException("ERROR: Res " + curPDBNum + " in rot/trans strand or BBFreeBlock but not flexible!");
            }

            resList.add(curRes);
        }

        return resList;
    }
    
    private ArrayList<DegreeOfFreedom> strandMotionDOFs(ArrayList<String[]> moveableStrands,
            ArrayList<String> flexibleRes){
        //Generate all the strand rotation/translation DOFs,
        //given the termini of the moving strands
        
        ArrayList<DegreeOfFreedom> ans = new ArrayList<>();
        
        for(String[] termini : moveableStrands){
            ArrayList<Residue> curStrandRes = resListFromTermini(termini, flexibleRes);
            MoveableStrand str = new MoveableStrand(curStrandRes);
            ans.addAll(str.getDOFs());
        }
        
        return ans;
    }
    
    
    private ArrayList<DegreeOfFreedom> strandDOFsForRes(Residue res, ArrayList<DegreeOfFreedom> strandDOFs){
        //List the strandDOFs that move res
        ArrayList<DegreeOfFreedom> ans = new ArrayList<>();
        
        for(DegreeOfFreedom dof : strandDOFs){
            MoveableStrand str;
            if(dof instanceof StrandRotation)
                str = ((StrandRotation)dof).getMoveableStrand();
            else//must be translation
                str = ((StrandTranslation)dof).getMoveableStrand();
            
            if(str.getResidues().contains(res))
                ans.add(dof);
        }
        
        return ans;
    }    
    
    //newPostList is an ArrayList that contains the positions to be merges 
    //(i.e., [1,3,4] means merge positions 1, 3, and 4
    public void mergePositions(ArrayList<Integer> posToMerge){
        //newPosList will be passed to helper function mergePosition
        //In this example it will be [[0],[1,3,4],[2],[5]]
        ArrayList<ArrayList<Integer>> newPosList = new ArrayList<>();
        
        //index into newPostList of new merged position
        //if newPosList is [[0],[1,3,4],[2],[5]], then newPosIndex = 1;
        int newPosIndex= -1;
        //new merge position consisting of multiple old positions
        ArrayList<Integer> newPos = new ArrayList<>();
        
        //iterate overall positions
        for (int i = 0; i<this.numPos; i++){
            //if this position is not going to be merged create a arraylist
            //that contains only this element
            if (! posToMerge.contains(i)){
                ArrayList<Integer> pos = new ArrayList<>();
                pos.add(i);
                newPosList.add(pos);
            }
            //if this position is going to be merged
            else{
                //check if newPosIndex has been set
                if (newPosIndex > -1){
                    //if it has been, then we add this position to the arraylist
                    //at the index of newPosList specified by newPosIndex
                    newPosList.get(newPosIndex).add(i);
                }
                else{
                    //if it has not been set then we set it and create the arraylist
                    //in this example i = 1 will reach this step.
                    newPosIndex = i;
                    ArrayList<Integer> mergePos = new ArrayList<>();
                    mergePos.add(i);
                    newPosList.add(mergePos);
                }
            }
        }
        mergePosition(newPosList);
    }
    
    //MergePositions is called by mergePositions to actually do the merging, 
    //Input: Arraylist of new position organization (i.e., [[0],[1,3],[2],[4]] 
    //       means merge positions 1,3
    //Creates a new PosFlex and overrides old PosFlex
    public void mergePosition(ArrayList<ArrayList<Integer>> newPosList){
        
        //Begin: Create new posFlex
        ArrayList<PositionConfSpaceSuper> newPosFlex = new ArrayList<>();
        int newPosFlexSize = newPosList.size();

        //Find original PosisitionConfSpaceSupers to merge
        for (int posIndex=0; posIndex<newPosFlexSize; posIndex++){
            ArrayList<PositionConfSpaceSuper> toMergePosConfSpaceSuper = new ArrayList<>();
            ArrayList<Integer> posList = newPosList.get(posIndex);
            for(int posConfSpaceSuperIndex : posList ){
                PositionConfSpaceSuper originalPosConfSpace = this.posFlex.get(posConfSpaceSuperIndex);
                toMergePosConfSpaceSuper.add(originalPosConfSpace);
            }
            //Merge PositionConfSpaces
            PositionConfSpaceSuper newPosConfSpace = mergePosConfs(toMergePosConfSpaceSuper);
            //Add new merged PositionConfSpace
            newPosFlex.add(newPosConfSpace);
        }
        //End: Create new posFlex
        
        this.posFlex = newPosFlex;
        this.numPos = newPosFlex.size();
    }

    //Given a list of PosConfSupers this method merges them all into one PosConfSuper        
    public PositionConfSpaceSuper mergePosConfs(ArrayList<PositionConfSpaceSuper> posConfSpaceList){
        PositionConfSpaceSuper newPosConfSpace;
        if(posConfSpaceList.size() == 1){//If size is 1 there is nothing to merge, we just return the position
            newPosConfSpace = posConfSpaceList.get(0);
        }
        else if(posConfSpaceList.size() == 2){//If only two, we can use the constructor
            newPosConfSpace = new PositionConfSpaceSuper(posConfSpaceList.get(0), posConfSpaceList.get(1));
        }
        else{//If more than two, we recurse
            ///This is not optimal in terms of speed, but for simplicity we will keep it now
            ArrayList<PositionConfSpaceSuper> toMerge = new ArrayList<>();//ArrayList of size 2 that we can then merge
            toMerge.add(new PositionConfSpaceSuper(posConfSpaceList.get(0), posConfSpaceList.get(1)));
            for (int i=2; i<posConfSpaceList.size(); i++){
                toMerge.add(posConfSpaceList.get(i));
            }
            newPosConfSpace = mergePosConfs(toMerge);
        }
        return newPosConfSpace;
    }
    
    private void standardizeMutatableRes(ArrayList<ArrayList<String>> allowedAAs, ArrayList<String> flexibleRes){
        //"mutate" all mutatable residues to the template version of their residue type
        //this ensures that the non-adjustable DOFs (bond angles, etc.) will be as in the template
        //(for consistency purposes)
        for(int pos=0; pos<numPos; pos++){
            Residue res = m.getResByPDBResNumber( flexibleRes.get(pos) );
            String resName = res.template.name;
            
            if(EnvironmentVars.resTemplates.getTemplateForMutation(resName, res, false) != null){
                //mutation to current residue type is possible, i.e., this position is mutatable
            //if(HardCodedResidueInfo.canMutateTo(res.template)){//this messes up for N-term mutations
                
                ResidueTypeDOF mutDOF = mutDOFs.get(pos);
                
                for(String allowedAA : allowedAAs.get(pos)){
                    if(allowedAA.equalsIgnoreCase("PRO")){
                        //mutating to and from PRO alters the sidechain a little, including idealization
                        //(because the sidechain connects in two places, it behaves a little different)
                        //Once we mutate to Pro and back we have a consistent conformation
                        mutDOF.mutateTo("PRO");
                        break;
                    }
                }                
                
                mutDOF.mutateTo(resName);
            }
        }
    }
  
    static ArrayList<DegreeOfFreedom> mutableResDOFs(Residue res, ArrayList<String> allowedAAs) {//mutation and dihedral confDOFs for the specified position

        ArrayList<DegreeOfFreedom> ans = new ArrayList<>();
        //assuming for now that this is a single-residue position...if multiple just add the following confDOFs for each residue

        //mutation DOF
        ans.add(new ResidueTypeDOF(res));

        int maxNumDihedrals = 0;//we need to create enough dihedral confDOFs for the allowed AA type with the most dihedrals

        for (String AAType : allowedAAs) {
            maxNumDihedrals = Math.max(maxNumDihedrals, EnvironmentVars.resTemplates.numDihedralsForResType(AAType));
        }

        for (int dih = 0; dih < maxNumDihedrals; dih++) {
            ans.add(new FreeDihedral(res, dih));
        }

        return ans;
    }
    
    //Get number of superRCs or RCs at each position
    public int[] getNumRCsAtPos(){
        int[] numAllowed = new int[numPos];
        
        for(int pos=0; pos<numPos; pos++){
            numAllowed[pos] = posFlex.get(pos).superRCs.size();
        }
        return numAllowed;
    }
    
    
     public double minimizeEnergy(int[] conf, EnergyFunction efunc, String outputPDBFile){
        //minimize the energy of a conformation, within the DOF bounds indicated by conf (a list of RCs)
        //return the minimized energy
        //if outputPDBFile isn't null, then output the minimized conformation to that file
        
        SuperRCTuple RCs = new SuperRCTuple(conf);
        MolecEObjFunction energy = new MolecEObjFunction(efunc,this,RCs);
        
        DoubleMatrix1D optDOFVals;
        
        if(energy.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
            Minimizer min = new CCDMinimizer(energy,false);
            //with the generic objective function interface we can easily include other minimizers though

            
            //DEBUG!!!!!  Timing pre-minimization without PB
            /*ArrayList<EnergyFunction> terms = ((MultiTermEnergyFunction)energy.getEfunc()).getTerms();
            ArrayList<Double> coeffs = ((MultiTermEnergyFunction)energy.getEfunc()).getCoeffs();
            if( terms.get(terms.size()-1) instanceof PoissonBoltzmannEnergy ){
                PoissonBoltzmannEnergy pbe = (PoissonBoltzmannEnergy)terms.remove(terms.size()-1);
                double pbcoeff = coeffs.remove(terms.size()-1);
                long startTime = System.currentTimeMillis();
                DoubleMatrix1D startDOFVals = min.minimize();
                long time1 = System.currentTimeMillis();
                terms.add(pbe);
                coeffs.add(pbcoeff);
                ((CCDMinimizer)min).setInitVals(startDOFVals);
                optDOFVals = min.minimize();
            }
            else//NON-DEBUG!*/
                optDOFVals = min.minimize();
        }
        else//molecule is already in the right, rigid conformation
            optDOFVals = DoubleFactory1D.dense.make(0);
        
        double minE = energy.getValue(optDOFVals);//this will put m into the minimized conformation
        
        if(outputPDBFile!=null)
            PDBFileWriter.writePDBFile(m, outputPDBFile, minE);
        
        return minE;
    }
}
