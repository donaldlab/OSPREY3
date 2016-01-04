/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.comets.LME;
import edu.duke.cs.osprey.astar.kadee.KaDEETree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.SearchProblemSuper;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.StringParsing;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 *
 * @author hunternisonoff
 */
public class KaDEEDoer {

    ConfigFileParser cfp;

    // A search problem for Super RCs.
    //0: Bound 
    //1: Mutable UnBound
    //2: Mutable Bound
    SearchProblem[] searchSpaces;
    SearchProblemSuper[] searchSpaceSupers;
    
    LME objFcn; //objective function for the KaDEE search
    LME[] constraints; //constraints for the KaDEE searchf
    int numTreeLevels; //number of mutable positions
    final int numStates = 2; //number of states considered
    ArrayList<ArrayList<String>> AATypeOptions = null; //AA types allowed at each mutable position
    ArrayList<ArrayList<Integer>> mutable2StatePosNums = new ArrayList<>();
    
    double Ew; // energy window for enumerating conformations: 0 for just GMEC
    double I0 = 0; // initial value of iMinDEE pruning interval
    boolean doIMinDEE;//do iMinDEE

    boolean useContFlex;
    //boolean enumByPairwiseLowerBound;//are we using a search method that 
    //enumerates by pairwise lower bound (minDEE-style)
    //or by (possibly approximated) true energy (rigid, EPIC, or tuple expander)?

    boolean outputGMECStruct;//write the GMEC structure to a PDB file

    boolean useEPIC = false;
    boolean useTupExp = false;

    boolean checkApproxE = true;//Use the actual energy function to evaluate
    //each enumerated conformation rather than just relying on the EPIC or tup-exp approximation

    boolean useEllipses = false;

    ExpFunction ef = new ExpFunction();
    double constRT = PoissonBoltzmannEnergy.constRT;

    public KaDEEDoer(ConfigFileParser cfp) {
        this.cfp = cfp;
        Ew = cfp.params.getDouble("Ew", 0);
        doIMinDEE = cfp.params.getBool("imindee", false);
        if (doIMinDEE) {
            I0 = cfp.params.getDouble("Ival", 5);
        }

        useContFlex = cfp.params.getBool("doMinimize", false);
        if (useContFlex || doIMinDEE) {
            throw new RuntimeException("Continuous Flexibility Not Yet Supported with KaDEE");
        }

        useTupExp = cfp.params.getBool("UseTupExp", false);
        useEPIC = cfp.params.getBool("UseEPIC", false);

        checkApproxE = cfp.params.getBool("CheckApproxE", true);

        if (doIMinDEE && !useContFlex) {
            throw new RuntimeException("ERROR: iMinDEE requires continuous flexibility");
        }

        outputGMECStruct = cfp.params.getBool("OUTPUTGMECSTRUCT", false);

        useEllipses = cfp.params.getBool("useEllipses", false);
    }

    /**
     * Run KaDEE and compute the sequence with the highest K* score. TODO: For
     * now this class only computes the partition function for the sequence with
     * the highest K* score.
     */
    void doKaDEE() {
        double curInterval = I0;//For iMinDEE.  curInterval will need to be an upper bound
        this.searchSpaces = cfp.getMSDSearchProblems();
        KaDEETree tree = setupKaDEETree();

        
        int[] seq1 = tree.nextConf();
        
        exhaustiveKaDEESearch();
    }

    //getGMEC from lower bounds
    private double calcGMEC(SearchProblem aSearchSpace) {
        SearchProblem searchSpace = aSearchSpace;
        ConfSearch search = new ConfTree(searchSpace);
        int[] conf = search.nextConf();
        double E = searchSpace.lowerBound(conf);
        return E;
    }

    //Loads energy matrices and prune for COMETS
    private void loadEMatandPruneComets(double pruningInterval) {
        cfp.params.setValue("TYPEDEP", "TRUE");
        for (int state = 0; state < searchSpaces.length; state++) {
            SearchProblem searchProblem = this.searchSpaces[state];

            System.out.println("Precomputing Energy Matrix for " + searchProblem.name + " state");
            searchProblem.loadEnergyMatrix();

            System.out.println("Initializing Pruning for " + searchProblem.name + " state");
            initializePruning(searchProblem);
            PruningControl pruning = cfp.setupPruning(searchProblem, pruningInterval, useEPIC, useTupExp);
            pruning.prune();
        }
    }

    //Given three search problems (Bound, UnBound Prot, Unbound Lig) this function
    //sets up the KaDEE tree.
    //The nonmutable unbound state is added and used just as a constant to the objective function
    private KaDEETree setupKaDEETree() {

        //For each state, for each position, this contains a list of allowed 
        //amino acids at that position
        ArrayList<ArrayList<ArrayList<String>>> allowedAAsPerState = new ArrayList<>();
        for (int state = 0; state < searchSpaces.length; state++) {
            allowedAAsPerState.add(getAllowedAA(state));
        }

        //Load the energy matrices and do pruning
        double pruningInterval = 0.0;
        loadEMatandPruneComets(pruningInterval);
        /*
         for (SearchProblemSuper searchProblem : searchSpaces){
         searchProblem.emat.setConstTerm(0.0);
         }
         */
        //For each state this arraylist gives the mutable pos nums of that state
        this.mutable2StatePosNums = handleMutable2StatePosNums(allowedAAsPerState);

        //determine which states are mutable and which are non-mutable
        boolean[] stateIsMutable = new boolean[this.searchSpaces.length];
        int numMutableState = 0;
        int numNonMutableState = 0;
        for (int state = 0; state < searchSpaces.length; state++) {
            stateIsMutable[state] = !(mutable2StatePosNums.get(state).isEmpty()); //If not empty, it is mutable
            if (stateIsMutable[state]) {
                numMutableState++;
            } else {
                numNonMutableState++;
            }
        }

        //Get Mutable States
        SearchProblem[] mutableStates = new SearchProblem[numMutableState];
        //For each MUTABLE state, this arraylist gives the mutable pos nums of that state
        //This is what will go into the COMETS tree
        ArrayList<ArrayList<Integer>> mutableState2StatePosNumList = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<String>>> mutableStateAllowedAAs = new ArrayList<>();
        int mutableStateIndex = 0;

        SearchProblem nonMutableState = searchSpaces[1];
        double unboundLigandGMECEnergy = 0.0;
        for (int state = 0; state < searchSpaces.length; state++) {
            if (stateIsMutable[state]) {
                mutableStates[mutableStateIndex] = searchSpaces[state];
                mutableState2StatePosNumList.add(mutable2StatePosNums.get(state));
                mutableStateAllowedAAs.add(allowedAAsPerState.get(state));
                mutableStateIndex++;
            } else {
                nonMutableState = searchSpaces[state];
                //For the const term of the LME objective function
                unboundLigandGMECEnergy = getGMECEnergyProtein(nonMutableState);
            }
        }
        //For const term of LME objective function
        int numStatesForCOMETS = mutableStates.length;
        this.numTreeLevels = getNumMutablePos(mutableState2StatePosNumList);
        this.AATypeOptions = handleAATypeOptions(mutableStateAllowedAAs);
        this.objFcn = new LME(new double[]{1, -1}, -unboundLigandGMECEnergy, 2);
        this.constraints = new LME[0];
        int numMaxMut = -1;
        String[] wtSeq = null;

        //Convert to ArrayList<>
        ArrayList<ArrayList<Integer>> mutableState2StatePosNum = new ArrayList<>();
        for (List<Integer> mutable2PosNum : mutableState2StatePosNumList) {
            ArrayList<Integer> converted = new ArrayList(mutable2PosNum);
            mutableState2StatePosNum.add(converted);
        }
        KaDEETree tree = new KaDEETree(numTreeLevels, objFcn, constraints, AATypeOptions, numMaxMut, wtSeq, mutableStateIndex, mutableStates, nonMutableState, mutableState2StatePosNum);
        return tree;
    }

    private ArrayList<ArrayList<String>> getAllowedAA(int state) {
        ArrayList<ArrayList<String>> complexAllowedAAs = cfp.getAllowedAAs();
        ArrayList<String> complexFlexRes = cfp.getFlexRes();
        int numFlexRes = complexFlexRes.size();
        int beginPos = -1;
        int endPos = -1;
        if (state == 0) {
            beginPos = 0;
            endPos = numFlexRes;
        } else if (state == 1) {
            beginPos = 0;
            endPos = this.searchSpaces[1].confSpace.posFlex.size();
        } else {
            beginPos = this.searchSpaces[1].confSpace.posFlex.size();
            endPos = beginPos + this.searchSpaces[2].confSpace.posFlex.size();
        }
        ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
        Molecule wtMolec = PDBFileReader.readPDBFile(cfp.params.getValue("PDBNAME"));

        for (int posNum = beginPos; posNum < endPos; posNum++) {
            ArrayList<String> currentAAOptions = complexAllowedAAs.get(posNum);
            if (cfp.params.getBool("AddWT", true)) {
                Residue res = wtMolec.getResByPDBResNumber(complexFlexRes.get(posNum));
                if (!StringParsing.containsIgnoreCase(complexFlexRes, res.template.name)) {
                    currentAAOptions.add(res.template.name);
                }
            }
            if (currentAAOptions.isEmpty()) {
                throw new RuntimeException("ERROR: No AAtype for Position: " + posNum);
            }
            allowedAAs.add(currentAAOptions);
        }
        return allowedAAs;
    }

    private String[] getWTSequence(List<Integer> complexMut2StatePosNum) {
        String[] wtSequence = new String[complexMut2StatePosNum.size()];
        ArrayList<String> complexFlexRes = cfp.getFlexRes();
        Molecule wtMolec = PDBFileReader.readPDBFile(cfp.params.getValue("PDBNAME"));

        for (int posIndex = 0; posIndex < complexMut2StatePosNum.size(); posIndex++) {
            Integer mutPosNum = complexMut2StatePosNum.get(posIndex);
            String mutPos = complexFlexRes.get(mutPosNum);
            wtSequence[posIndex] = wtMolec.getResByPDBResNumber(mutPos).template.name;
        }
        return wtSequence;
    }

    //for each state, returns a list of mutable pos numbers. 
    //If the state has no mutable positions an empty list is returned
    //Thus, we can use this to decide which states go in the comets search
    private ArrayList<ArrayList<Integer>> handleMutable2StatePosNums(ArrayList<ArrayList<ArrayList<String>>> allowedAAPerState) {
        ArrayList<ArrayList<Integer>> mutable2StatePosNum = new ArrayList<>();

        for (int state = 0; state < allowedAAPerState.size(); state++) {
            ArrayList<ArrayList<String>> allowedAAForState = allowedAAPerState.get(state);
            ArrayList<Integer> mutablePositionsForState = getMutablePosNums(allowedAAForState);
            mutable2StatePosNum.add(mutablePositionsForState);
        }

        return mutable2StatePosNum;
    }

    //Given the allowed AAs at a state, get the mutable position numbers
    private ArrayList<Integer> getMutablePosNums(ArrayList<ArrayList<String>> allowedAAForState) {
        ArrayList<Integer> mutablePosNum = new ArrayList<>();
        for (int posNum = 0; posNum < allowedAAForState.size(); posNum++) {
            if (allowedAAForState.get(posNum).size() > 1) {
                mutablePosNum.add(posNum);
            }
        }
        return mutablePosNum;
    }

    private int getNumMutablePos(ArrayList<ArrayList<Integer>> mutable2StatePosNum) {
        int numMutable = -1;
        for (int state = 0; state < mutable2StatePosNum.size(); state++) {
            int numMutableAtState = mutable2StatePosNum.get(state).size();
            numMutable = Math.max(numMutable, numMutableAtState);
        }
        return numMutable;
    }

    //Make sure each mutable state has the same Allowed AA
    //Return Allowed AA for all mutable positions
    private ArrayList<ArrayList<String>> handleAATypeOptions(ArrayList<ArrayList<ArrayList<String>>> mutableStateAllowedAAs) {
        ArrayList<ArrayList<String>> AATypeOptions = mutableStateAllowedAAs.get(0).stream().filter(aaTypes -> aaTypes.size() > 1).collect(Collectors.toCollection(ArrayList::new));

        for (int state = 1; state < mutableStateAllowedAAs.size(); state++) {
            ArrayList<ArrayList<String>> AATypesForState = mutableStateAllowedAAs.get(state).stream().filter(aaTypes -> aaTypes.size() > 1).collect(Collectors.toCollection(ArrayList::new));
            if (AATypesForState.size() != AATypeOptions.size()) {
                throw new RuntimeException("ERROR: Different Number of Mutable Positions between Bound and Unbound");
            }
            for (int posNum = 0; posNum < AATypesForState.size(); posNum++) {
                ArrayList<String> AATypesForPos = AATypesForState.get(posNum);
                for (int aaIndex = 0; aaIndex < AATypesForPos.size(); aaIndex++) {
                    if (!(AATypeOptions.get(posNum).contains(AATypesForPos.get(aaIndex)))) {
                        throw new RuntimeException("ERROR: AAType Different for Bound and Unbound Mutable Residues");
                    }
                }
            }
        }
        return AATypeOptions;
    }

    private double getGMECEnergyProtein(SearchProblem searchProblem) {
        ConfTree confTree = new ConfTree(searchProblem);
        int[] conf = confTree.nextConf();
        RCTuple rcTup = new RCTuple(conf);
        double unboundLigandE = searchProblem.emat.getInternalEnergy(rcTup) + searchProblem.emat.getConstTerm();
        return unboundLigandE;
    }

    private void initializePruning(SearchProblem searchProblem) {
        //Creates an efficient competitor pruning matrix
        searchProblem.competitorPruneMat = null;
        System.out.println("PRECOMPUTING COMPETITOR PRUNING MATRIX");
        //prune with 0 interval, anything that survives will be added as a competitor
        PruningControl compPruning = cfp.setupPruning(searchProblem, 0, false, false);
        compPruning.setOnlyGoldstein(true);
        compPruning.prune();
        searchProblem.competitorPruneMat = searchProblem.pruneMat;
        searchProblem.pruneMat = null;
        System.out.println("COMPETITOR PRUNING DONE");
    }

    //For quality control, it's good to be able to check KaDEE results by exhaustive search
    void exhaustiveKaDEESearch() {

        System.out.println();
        System.out.println("CHECKING KaDEE RESULT BY EXHAUSTIVE SEARCH");
        System.out.println();

        ArrayList<String[]> seqList = listAllSeqs();
        int numSeqs = seqList.size();
        double stateGMECs[][] = new double[numSeqs][3];

        for (int state = 0; state < this.numStates; state++) {
            for (int seqNum = 0; seqNum < numSeqs; seqNum++) {
                String[] sequence = seqList.get(seqNum);
                SearchProblem searchProb = (SearchProblem) ObjectIO.deepCopy(this.searchSpaces[state]);
                handleStatePruning(searchProb, sequence, this.mutable2StatePosNums.get(state));
                stateGMECs[seqNum][state] = getMAP(searchProb);
            }
        }

        //now find the best sequence and obj fcn value
        int topSeqNum = -1;
        double bestSeqScore = Double.POSITIVE_INFINITY;

        for (int seqNum = 0; seqNum < numSeqs; seqNum++) {
            boolean constrSatisfied = true;
            for (LME constr : constraints) {
                if (constr.eval(stateGMECs[seqNum]) > 0) {
                    constrSatisfied = false;
                }
            }

            if (constrSatisfied) {
                double seqScore = objFcn.eval(stateGMECs[seqNum]);
                if (seqScore < bestSeqScore) {
                    bestSeqScore = seqScore;
                    topSeqNum = seqNum;
                }
            }
        }

        System.out.println();
        if (topSeqNum == -1) {
            System.out.println("EXHAUSTIVE SEARCH FINDS NO CONSTR-SATISFYING SEQUENCES");
        } else {
            System.out.println("EXHAUSTIVE MULTISTATE BEST SCORE: " + bestSeqScore + " SEQUENCE: ");
            for (String aa : seqList.get(topSeqNum)) {
                System.out.println(aa);
            }
        }
        System.out.println();
    }

    private void handleStatePruning(SearchProblem searchProb, String[] sequence, ArrayList<Integer> mutable2StatePosNum){
        if (mutable2StatePosNum.size() != sequence.length){
            throw new RuntimeException("ERROR: Length of Sequence Not Equal to NumMutablePositions. Cannot complete exhaustive search");
        }
        
        for (int i=0; i<mutable2StatePosNum.size(); i++){
            int posNum = mutable2StatePosNum.get(i);
            String AAtype = sequence[i];
            for (int rc : searchProb.pruneMat.unprunedRCsAtPos(posNum)){
                String rcAAType = searchProb.confSpace.posFlex.get(posNum).RCs.get(rc).AAType;
                
                if (!rcAAType.equalsIgnoreCase(AAtype)){
                    searchProb.pruneMat.markAsPruned(new RCTuple(posNum,rc));
                }
            }
        }
        
    }
    
    private ArrayList<String[]> listAllSeqs() {
        //List all possible sequence for the mutable residues,
        //based on AATypeOptions
        return listAllSeqsHelper(0);
    }

    private ArrayList<String[]> listAllSeqsHelper(int mutPos) {
        //List all partial sequences for the subset of mutable positions
        //starting at mutPos and going to the last mutable position
        ArrayList<String[]> ans = new ArrayList<>();

        if (mutPos == numTreeLevels) {
            ans.add(new String[0]);
        } else {
            ArrayList<String[]> subList = listAllSeqsHelper(mutPos + 1);
            for (String AAType : AATypeOptions.get(mutPos)) {
                for (String[] subSeq : subList) {
                    String seq[] = new String[numTreeLevels - mutPos];
                    System.arraycopy(subSeq, 0, seq, 1, numTreeLevels - mutPos - 1);
                    seq[0] = AAType;
                    ans.add(seq);
                }
            }
        }

        return ans;
    }


    private double getMAP(SearchProblem searchSpace) {
        ConfTree confTree = new ConfTree(searchSpace);

        if (searchSpace.contSCFlex) {
            throw new RuntimeException("Continuous Flexibility Not Yet Supported in KaDEE");
        }

        int[] MAPconfig = confTree.nextConf();
        double E = searchSpace.emat.getInternalEnergy(new RCTuple(MAPconfig));
        return E;
    }
}
