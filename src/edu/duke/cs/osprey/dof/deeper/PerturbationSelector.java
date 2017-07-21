/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof.deeper;

import java.util.ArrayList;
import java.util.TreeSet;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.dof.deeper.perts.Perturbation;
import edu.duke.cs.osprey.multistatekstar.ResidueTermini;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residue.SecondaryStructure;

/**
 *
 * @author mhall44
 */
public class PerturbationSelector {
    
    
    //parameters needed
    String startingPertFile; 
    boolean onlyStarting;
    double maxShearParam, maxBackrubParam;//upper limit on these (by default) single-interval
    //perturbation parameters.  We also allow this value in the negative direction
    boolean selectLCAs;
    boolean doRamaCheck;//Check if proposed perturbation states are Ramachandran allowed
    ArrayList<String> flexibleRes;//PDB numbers of flexible residues
    
    Strand strand;//generated from provided PDB file
    
    
    
    
    //Stuff generated and used in the selection process
    ArrayList<Perturbation> perts;//the perturbations, implemented in m
    PertSet ps;//the perturbation set we're generating
    
    ArrayList<TreeSet<String>> resMovedByPert;//for each perturbation, the set of residues
    //whose conformations depend on its parameter value
    
    ArrayList<double[]> defaultLCAIntervals, defaultShearIntervals, defaultBackrubIntervals;
    
    int failingPertIndex = -1;//when we are trying combinations of perturbations for a residue,
    //a value other than -1 means perturbation # failingPertIndex could not be applied
    //so we should avoid all other combination of intervals that start the same way

    
    
    public PerturbationSelector(String startingPertFile, boolean onlyStarting, 
            double maxShearParam, double maxBackrubParam, boolean selectLCAs, 
            ArrayList<String> flexibleRes, String PDBFile, ResidueTermini termini, 
            boolean doRamaCheck, ResidueTemplateLibrary templateLib) {
        
        this.startingPertFile = startingPertFile;
        this.onlyStarting = onlyStarting;
        this.maxShearParam = maxShearParam;
        this.maxBackrubParam = maxBackrubParam;
        this.selectLCAs = selectLCAs;
        this.flexibleRes = flexibleRes;
        
        strand = new Strand.Builder(PDBIO.readFile(PDBFile))
            .setTemplateLibrary(templateLib)
            .setResidues(termini)
            .build();
    }
    
    
    
    public PertSet selectPerturbations(ResidueTermini termini){
        ps = new PertSet();
        
        if( !startingPertFile.equalsIgnoreCase("none") ){//lists perturbations to include
            //that would not be automatically selected.  (Probably partial structure switches)
            if( !ps.loadPertFile(startingPertFile,false, termini) ){
                throw new RuntimeException( "ERROR: Can't find starting perturbation file "
                        + startingPertFile );
            }
        }
        
        if(!onlyStarting){//select other perturbations as appropriate
            //giving them default intervals
            autogeneratePerturbations();
        }
        
        calcResMovedByPert();
        
        mutateFlexResToGly();

        perts = ps.makePerturbations(strand.mol);
        ps.pertStates = new ArrayList<>();
        
        //OK now figure out which states are available for each residue
        for(int pos=0; pos<flexibleRes.size(); pos++){
            ArrayList<ArrayList<int[]>> resPertStates = new ArrayList<>();
            
            ArrayList<Integer> pertIndices = pertIndicesForPos(pos);
            ArrayList<int[]> state = new ArrayList<>();
            for(int pertInd : pertIndices)
                state.add(new int[]{pertInd,0});//start in first interval for each pert
            
            while(state!=null){//go through all combinations of intervals
                if(isStateReasonable(state,pos)){
                    resPertStates.add(state);
                }
                
                state = nextPossibleState(state);
                //increment to next state, or skip over clearly impossible ones
                //(based on failingPertIndex)
            }
            
            ps.pertStates.add(resPertStates);
        }
        
        removeIncompatiblePertStates(ps.pertStates);
        
        return ps;
    }
    
    
    void mutateFlexResToGly(){
        for(String resNum : flexibleRes){
       
        	// AAO 2016: mutation assumes residue is an amino acid. throws an exception otherwise
            Residue res = strand.mol.getResByPDBResNumber(resNum);
            if(HardCodedResidueInfo.hasAminoAcidBB(res) && !res.fullName.startsWith("FOL"))
            	new ResidueTypeDOF(strand.templateLib, res).mutateTo("GLY");
        }
    }
    
    
    private void calcResMovedByPert(){
        
        resMovedByPert = new ArrayList<>();
        
        for(int pertIndex=0; pertIndex<ps.pertTypes.size(); pertIndex++){
            
            TreeSet<String> curPertRes = new TreeSet<>();

            curPertRes.addAll(ps.resNums.get(pertIndex));//residue directly affected by pertIndex
            
            for(int pertIndex2=pertIndex+1; pertIndex2<ps.pertTypes.size(); pertIndex2++){
                //for any perturbation coming after pertIndex,
                //if our perturbation (pertIndex) can change its starting conformation,
                //then it can change the conformations of all residues it affects directly
                ArrayList<String> pert2Res = ps.resNums.get(pertIndex2);
                for(String resNum : pert2Res){
                    if(curPertRes.contains(resNum)){
                        curPertRes.addAll(pert2Res);
                        break;
                    }
                }
            }
            
            resMovedByPert.add(curPertRes);
        }
    }
    
    
    ArrayList<Integer> pertIndicesForPos(int pos){
        //which perturbations (indexed in perts) have an effect on the conformation of flexible position pos?
        ArrayList<Integer> ans = new ArrayList<>();
        String resNum = flexibleRes.get(pos);
        
        for(int pertIndex=0; pertIndex<perts.size(); pertIndex++){
            if( resMovedByPert.get(pertIndex).contains(resNum) )
                ans.add(pertIndex);
        }
        
        return ans;
    }

    
    private void removeIncompatiblePertStates(ArrayList<ArrayList<ArrayList<int[]>>> pertStates){
        //In order for a residue perturbation state to be possible,
        //it must be compatible with at least one perturbation state at each other residue
        //Remove those that aren't
        
        int numPos = pertStates.size();
        boolean prunedStates[][] = getIncompatiblePertStates(pertStates);
        
        for(int pos=0; pos<numPos; pos++){
            //remove walking backwards to avoid messing up indexing
            int resNumStates = pertStates.get(pos).size();
            
            for(int state=resNumStates-1; state>=0; state--){
                if(prunedStates[pos][state])
                    pertStates.get(pos).remove(state);
            }
        }
    }
    
    
    private boolean[][] getIncompatiblePertStates(ArrayList<ArrayList<ArrayList<int[]>>> pertStates){
        //What states are pruned (based on incompatibility) at each residue position?
        
        int numPos = pertStates.size();
        
        boolean prunedStates[][] = new boolean[numPos][];//Which states are pruned
        //at each residue position
        for(int pos=0; pos<numPos; pos++)
            prunedStates[pos] = new boolean[pertStates.get(pos).size()];

        
        boolean done = false;

        while( !done ){//We iterate until no more states can be removed

            done = true;

            for (int curPos=0; curPos<numPos; curPos++){

                int resNumStates = pertStates.get(curPos).size();

                for(int curState=0; curState<resNumStates; curState++){

                    if( ! prunedStates[curPos][curState] ){

                        ArrayList<int[]> state1 = pertStates.get(curPos).get(curState);

                        for (int altPos=0; altPos<numPos; altPos++){

                            if( altPos != curPos ){

                                boolean prune = true;//Indicates no perturbation state compatible with curRC has been found yet at altPos
                                int altResNumStates = pertStates.get(altPos).size();

                                for(int altState=0; altState<altResNumStates; altState++){
                                    
                                    if( ! prunedStates[altPos][altState] ){
                                    
                                        ArrayList<int[]> state2 = pertStates.get(altPos).get(altState);

                                        if( ! arePertStatesIncompatible(state1, state2) ){
                                            prune = false;
                                            break;
                                        }
                                    }
                                }

                                prunedStates[curPos][curState] = prunedStates[curPos][curState] || prune;
                            }
                        }

                        if( prunedStates[curPos][curState] )//Iterate again if anything was pruned
                            done = false;
                    }
                }
            }

        }

        return prunedStates;
    }
    
    
    
    public boolean arePertStatesIncompatible(ArrayList<int[]> state1, ArrayList<int[]> state2){
        //Two perturbation states at different residues are incompatible
        //if they have different intervals for the same perturbation
                
        for(int[] p1 : state1){
            for(int[] p2 : state2){
                if(p1[0]==p2[0]){//same perturbation
                    if(p1[1]!=p2[1])//different states
                        return true;
                }
            }
        }
        
        return false;//No incompatibility found, so we're good
    }
    
    
    private ArrayList<int[]> nextPossibleState(ArrayList<int[]> state){
        //Increment a residue perturbation state to the next possible combination
        //of intervals for its perturbations
        //Skip states if indicated by failingPertIndex
        //Return null if this is the last possible state
        
        if(state.isEmpty())//no perturbations: only one state possible (unperturbed)
            return null;
        
        //first, copy state
        ArrayList<int[]> ans = new ArrayList<>();
        for(int[] p : state)
            ans.add(p.clone());
        
        if(failingPertIndex==-1)//increment last perturbation's interval
            return incrementPertState(ans,ans.size()-1);
        else//increment failing perturbation (keeping current value would surely fail)
            return incrementPertState(ans,failingPertIndex);
    }
    
    
    private ArrayList<int[]> incrementPertState(ArrayList<int[]> state, int pertIndex){
        //Increment the interval in state of the given perturbation (indexed in state)
        //If maxed out, return null
        
        int curInterval = state.get(pertIndex)[1];//index (in pertIntervals) of current interval
        //allowed to range from 0 to maxInterval
        int maxInterval = ps.pertIntervals.get(state.get(pertIndex)[0]).size()-1;
        
        if(curInterval<maxInterval){//can increment here
            state.get(pertIndex)[1]++;
            return state;
        }
        else {//cannot...increment previous perturbation
            
            if(pertIndex==0)//no previous one...we've reached the last state
                return null;
            
            //revert all later perturbations to 0
            for(int index2=pertIndex; index2<state.size(); index2++)
                state.get(pertIndex)[1] = 0;
            
            return incrementPertState(state,pertIndex-1);
        }
    }
    
    
    private void initDefaultIntervals(){
        //initialize the default intervals for LCAs, shears, and backrubs
        
        defaultLCAIntervals = new ArrayList<>();
        //LCAs are discrete but can have up to 16 states
        for(double param=0; param<16; param++)
            defaultLCAIntervals.add(new double[] {param,param});
        
        defaultShearIntervals = new ArrayList<>();
        defaultShearIntervals.add(new double[] {-maxShearParam,maxShearParam});
        
        defaultBackrubIntervals = new ArrayList<>();
        defaultBackrubIntervals.add(new double[] {-maxBackrubParam,maxBackrubParam});
    }
    
    
    private void autogeneratePerturbations(){
        //generate the perturbations we want in ps, based on
        //stretches of consecutive flexible residues that could accommodate them
        
        initDefaultIntervals();

        //get lists of consecutive residue numbers that can be used to
        //make perturbations
        ArrayList<ArrayList<String>> consecTriplesBR = consecutiveFlexibleRes(3,"BACKRUB");
        ArrayList<ArrayList<String>> consecTriplesLCA = consecutiveFlexibleRes(3,"LOOP CLOSURE ADJUSTMENT");
        ArrayList<ArrayList<String>> consecQuads = consecutiveFlexibleRes(4,"SHEAR");

        if(selectLCAs){
            //start with these because they're discrete but make lots of states
            for(ArrayList<String> lcaRes : consecTriplesLCA){
                ps.pertTypes.add("LOOP CLOSURE ADJUSTMENT");
                ps.resNums.add(lcaRes);
                ps.pertIntervals.add(defaultLCAIntervals);
                ps.additionalInfo.add(null);//LCAs don't need additional info
            }
        }


        //ok now shears and backrubs
        for(ArrayList<String> shearRes : consecQuads){
            ps.pertTypes.add("SHEAR");
            ps.resNums.add(shearRes);
            ps.pertIntervals.add(defaultShearIntervals);
            ps.additionalInfo.add(null);
        }

        for(ArrayList<String> brRes : consecTriplesBR){
            ps.pertTypes.add("BACKRUB");
            ps.resNums.add(brRes);
            ps.pertIntervals.add(defaultBackrubIntervals);
            ps.additionalInfo.add(null);
        }

        //partial structure switches must be put in by startingPertFile
    }
    
    
    ArrayList<ArrayList<String>> consecutiveFlexibleRes(int resCount, String pertType){
        //Return subsets of the flexibleRes that are resCount consecutive residues
        //we'll also screen them based on the type of perturbation we want to put here
        
        ArrayList<Integer> startPos = new ArrayList<>();
        //which flexible residue positions are good for starting this perturbation?
        
        //Look for any sets of resCount residues in a row in flexMolResNum
        for( int pos=0; pos<flexibleRes.size() - resCount + 1; pos++ ){
            
            //OK so for pos to be a valid position to start a perturbation at,
            //it must the start of resCount flexible residues peptide-bonded
            //to each other, and they must have the right secondary structure
            
            //check bonding
            boolean correctBonding = true;
            
            for(int offset=1; offset<resCount; offset++){
                Residue res1 = strand.mol.getResByPDBResNumber(flexibleRes.get(pos+offset-1));
                Residue res2 = strand.mol.getResByPDBResNumber(flexibleRes.get(pos+offset));
                
                int Cindex = res1.getAtomIndexByName("C");
                int Nindex = res2.getAtomIndexByName("N");
                
                if(Cindex==-1 || Nindex==-1){
                    correctBonding = false;
                    break;
                }
                
                Atom C = res1.atoms.get(Cindex);
                Atom N = res2.atoms.get(Nindex);
                if(!C.bonds.contains(N)){
                    correctBonding = false;
                    break;
                }
            }
            
            if(correctBonding){
                //we have a stretch of bonded residues.  Check secondary structure.
                ArrayList<Residue> candidateRes = new ArrayList<>();
                for(int offset=0; offset<resCount; offset++){
                    candidateRes.add(strand.mol.getResByPDBResNumber(flexibleRes.get(pos+offset)));
                }
                if( secondaryStructureCorrect(candidateRes,pertType) ){
                    startPos.add(pos);
                }
            }
        }
        
        ArrayList<ArrayList<String>> ans = new ArrayList<>();

        for(int shift=0; shift<resCount; shift++){
        //We add perturbations in order of their starting molecule residue number modulo resCount,
        //then in order of the actual molecule residue number, to minimize extra RCs from overlaps
            for(int pos : startPos){
                if( pos % resCount == shift ){
                    
                    ArrayList<String> pertResNums = new ArrayList<>();
                    for(int offset=0; offset<resCount; offset++)
                        pertResNums.add(flexibleRes.get(pos+offset));

                    ans.add(pertResNums);
                }
            }
        }
        
        return ans;
    }
    
    
    private boolean secondaryStructureCorrect(ArrayList<Residue> resList, String pertType){
        //Do the residues have the right secondary structure to undergo this perturbation?
        //we'll require at least 3 helix residues for a shear; all 3 non-helix for backrub;
        //at least 2 loop for loop closure adjustment
        
        boolean good = true;
        
        if(pertType.equalsIgnoreCase("SHEAR")){
            int count = 0;
            for(Residue res : resList){
                if( res.secondaryStruct == SecondaryStructure.HELIX )
                    count++;
            }
            if(count < 3)
                good = false;
        }

        else if(pertType.equalsIgnoreCase("BACKRUB")){
            for(Residue res : resList){
                if( res.secondaryStruct == SecondaryStructure.HELIX )
                    good = false;
            }
        }

        else if(pertType.equalsIgnoreCase("LOOP CLOSURE ADJUSTMENT")){
            int count = 0;
            for(Residue res : resList){
                if( res.secondaryStruct == SecondaryStructure.LOOP )
                    count++;
            }
            if(count < 2)
                good = false;
        }
        
        else
            throw new RuntimeException("ERROR: Unrecognized perturbation: "+pertType);
        
        return good;
    }
    
    
    private boolean isStateReasonable(ArrayList<int[]> pertState, int pos){
        //Try to apply the perturbation state to see if it's possible,
        //and if so if it's OK Ramachandran-wise
        
        //if pertState all unperturbed then keep it
        boolean allUnperturbed = true;
        for(int[] p : pertState){
            if(p[1]!=0)
                allUnperturbed = false;
        }
        if(allUnperturbed){
            failingPertIndex = -1;
            return true;
        }
        
        
        //if needed build a structure based on first several perts
        //that will help us quickly rule out perts...
        
        double[][] backupCoords = backupFlexResCoords();//backup coordinates for all flexible residues
        
        for(int pertInd=0; pertInd<pertState.size(); pertInd++){
            int[] p = pertState.get(pertInd);
            double interval[] = ps.pertIntervals.get(p[0]).get(p[1]);
            double midVal = 0.5*(interval[0]+interval[1]);//we'll check pert feasibility at middle of interval
            if(midVal!=0){//there is a motion to perform
                if(!perts.get(p[0]).doPerturbationMotion(midVal)){//not needed for Pro flip??  or midval=0?
                    failingPertIndex = pertInd;
                    return false;
                }
            }
        }
        
        //OK if we get here the pert state is geometrically possible
        //Check Ramachandran for this residue...
        failingPertIndex = -1;//did not fail to apply perturbation
        boolean ok = ramaCheck(strand.mol.getResByPDBResNumber(flexibleRes.get(pos)));//just check res,
        //see if works for any Ramachandran category except Gly
        
        restoreFlexResCoords(backupCoords);
        return ok;
    }
    
    
    private double[][] backupFlexResCoords(){
        //coolect coordinates from all the flexible residues
        int numPos = flexibleRes.size();
        double[][] backup = new double[numPos][];
        for(int flexRes=0; flexRes<numPos; flexRes++){
            Residue res = strand.mol.getResByPDBResNumber(flexibleRes.get(flexRes));
            backup[flexRes] = res.coords.clone();
        }
        
        return backup;
    }
    
    
    private void restoreFlexResCoords(double[][] backup){
        //restore from backup
        int numPos = flexibleRes.size();
        for(int flexRes=0; flexRes<numPos; flexRes++){
            Residue res = strand.mol.getResByPDBResNumber(flexibleRes.get(flexRes));
            res.coords = backup[flexRes];
        }
    }
    
    
    boolean ramaCheck(Residue res){
        //We'll do Ramachandran check based on gly, since the residue could mutate to gly maybe
        //VDW energies will indicate which gly-allowed states are not cool for a given mutation
        if(!doRamaCheck)
            return true;
        
        boolean[] allowed = RamachandranChecker.getInstance().checkByAAType(res);
        return allowed[0];//denotes Gly
    }
}
