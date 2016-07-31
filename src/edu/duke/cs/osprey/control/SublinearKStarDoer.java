/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.comets.LME;
import edu.duke.cs.osprey.astar.partfunc.PartFuncTree;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.astar.seqkstar.SublinearKStarTree;
import edu.duke.cs.osprey.partitionfunctionbounds.TRBP;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.StringParsing;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 *
 * @author hmn5
 */
public class SublinearKStarDoer {

    ConfigFileParser cfp;

    // A search problem for Super RCs.
    //0: Bound 
    //1: Mutable UnBound
    //2: Mutable Bound
    SearchProblem[] searchSpaces;
    SearchProblem[] mutableSearchSpace;

    SearchProblem boundSearchSpace; //the search space of the protein-ligand bound complex
    SearchProblem unboundMutSearchSpace; //the search space of the unbound mutable (protein) 
    SearchProblem unboundNonMutSearchSpace; // the search space of the unbound non-mutable (ligand)

    LME objFcn; //objective function 
    LME[] gmecConstraints; //constraints on the GMEC 
    LME[] ensembleConstraints; // constraints on the free-energy/change in free-energy

    int numTreeLevels; //number of mutable positions

    final int numStates = 2; //number of states considered
    ArrayList<ArrayList<String>> AATypeOptions = null; //AA types allowed at each mutable position
    ArrayList<ArrayList<Integer>> mutable2StatePosNums = new ArrayList<>();

    boolean doExhaustive = false;
    boolean computeNonMutPF = true; //Compute the non-mutable partition function (for obj func const term)

    ExpFunction ef = new ExpFunction();
    double constRT = PoissonBoltzmannEnergy.constRT;

    public String[] bestSequence;

    public SublinearKStarDoer(ConfigFileParser cfp) {
        this.cfp = cfp;

        boolean useContinuousFlex = cfp.params.getBool("doMinimize");
        if (useContinuousFlex) {
            throw new RuntimeException("Continuous Flexibility Not Yet Supported");
        }

        boolean useTupExp = cfp.params.getBool("UseTupExp");
        boolean useEPIC = cfp.params.getBool("UseEPIC");
        if (useTupExp || useEPIC) {
            throw new RuntimeException("LUTE and EPIC Not Yet Supported");
        }
    }

    /**
     * Run SublinearKStar and compute the sequence with the highest K* score.
     * TODO: For now this class only computes the partition function for the
     * sequence with the highest K* score.
     * @param doExhaustive Should we do exhaustive search (for debugging/testing)
     */
    public void doSublinearKStar(boolean doExhaustive) {
        if (doExhaustive) {
            setupSublinearKStarTree();
            exhaustiveSublinearKStarSearch();
        } else {
            SublinearKStarTree tree = setupSublinearKStarTree();
            long startTime = System.currentTimeMillis();
            tree.nextConf();
            long totalTime = System.currentTimeMillis() - startTime;
            System.out.println("Sublinear KStar Finished in " + totalTime + " ms");
        }
    }

    public SublinearKStarTree setupSublinearKStarTree() {
        this.searchSpaces = cfp.getMSDSearchProblems();

        // Get the Bound Search Problem
        this.boundSearchSpace = getBoundSearchProblem(searchSpaces);
        // Get the amino acids allowed at each position for the bound search problem
        ArrayList<ArrayList<String>> allowedAAsPerBoundPos = getAllowedAA(this.boundSearchSpace);
        // Get the position numbers that have more than one amino acid for this search problem
        ArrayList<Integer> mutable2BoundStatePosNums = IntStream.range(0, allowedAAsPerBoundPos.size())
                .filter(pos -> allowedAAsPerBoundPos.get(pos).size() > 1)// filter pos numbers for more than one amino acid
                .boxed().collect(Collectors.toCollection(ArrayList::new));

        // Get the Unbound Mutable Search Problem
        this.unboundMutSearchSpace = getUnboundMutableSearchProblem(searchSpaces);
        // Get the amino acids allowed at each position for the unbound mutable search problem
        ArrayList<ArrayList<String>> allowedAsPerUnboundPos = getAllowedAA(this.unboundMutSearchSpace);
        // Get the position numbers that have more than one amino acid for this search problem
        ArrayList<Integer> mutable2UnboundStatePosNums = IntStream.range(0, allowedAsPerUnboundPos.size())
                .filter(pos -> allowedAAsPerBoundPos.get(pos).size() > 1)// filter pos numbers for more than one amino acid
                .boxed().collect(Collectors.toCollection(ArrayList::new));
        
        // Get the Unbound Non-Mutable Search Problem
        this.unboundNonMutSearchSpace = getUnboundNonMutableSearchProblem(searchSpaces);

        // Load the energy matrices and do pruning
        double pruningInterval = Double.POSITIVE_INFINITY;
        loadEMatandPrune(pruningInterval);
        
        // The constant term for the objective function is the unbound non-mutable partition function
        double unboundNonMutLogPartFunction;
        if (computeNonMutPF) {
            PartFuncTree pft = new PartFuncTree(this.unboundNonMutSearchSpace);
            unboundNonMutLogPartFunction = pft.computeEpsilonApprox(0.1);
        } else {
            unboundNonMutLogPartFunction = 0.0;
        }
        System.out.println("Unbound Non Mutable Log PF: " + unboundNonMutLogPartFunction);

        this.numTreeLevels = mutable2BoundStatePosNums.size();
        this.AATypeOptions = handleAATypeOptions(allowedAAsPerBoundPos, allowedAsPerUnboundPos);
        this.objFcn = new LME(new double[]{1, -1}, unboundNonMutLogPartFunction, 2);
        this.gmecConstraints = new LME[0];
        int numMaxMut = -1;
        String[] wtSeq = getWTSequence(mutable2BoundStatePosNums); 

        SublinearKStarTree tree = new SublinearKStarTree(numTreeLevels, objFcn, gmecConstraints, AATypeOptions, 
                numMaxMut, wtSeq, this.boundSearchSpace, this.unboundMutSearchSpace, this.unboundNonMutSearchSpace, 
                mutable2BoundStatePosNums, mutable2UnboundStatePosNums);

        return tree;
    }

    //Return the bound search problem
    private SearchProblem getBoundSearchProblem(SearchProblem[] searchProbs) {
        List<SearchProblem> searchProblems = Arrays.asList(searchProbs);
        //Sort in reverse order by number of residues
        SearchProblem boundSP = searchProblems.stream()
                .sorted((sp1, sp2) -> Integer.compare(sp2.confSpace.numPos, sp1.confSpace.numPos))
                .findFirst().get();
        return boundSP;
    }

    private SearchProblem getUnboundMutableSearchProblem(SearchProblem[] searchProbs) {
        List<SearchProblem> searchProblems = Arrays.asList(searchProbs);
        SearchProblem boundSP = getBoundSearchProblem(searchProbs);

        // get unbound search problems
        List<SearchProblem> unBoundSearchProblems = searchProblems.stream()
                .filter(searchProb -> !searchProb.name.equals(boundSP.name))
                .collect(Collectors.toList());

        // for each unbound search problem filter if all positions have only one amino acid
        List<SearchProblem> mutableUnboundSp = unBoundSearchProblems.stream()
                .filter(sp -> getAllowedAA(sp).stream().anyMatch(pos -> pos.size() > 1))
                .collect(Collectors.toList());

        // make sure we only have one mutable strand
        if (mutableUnboundSp.size() > 1) {
            throw new RuntimeException("ERROR: Only One Strand Can Currently Be Mutable");
        }
        return mutableUnboundSp.get(0);
    }

    private SearchProblem getUnboundNonMutableSearchProblem(SearchProblem[] searchProbs) {
        List<SearchProblem> searchProblems = Arrays.asList(searchProbs);
        SearchProblem boundSP = getBoundSearchProblem(searchProbs);

        // get unbound search problems
        List<SearchProblem> unBoundSearchProblems = searchProblems.stream()
                .filter(searchProb -> !searchProb.name.equals(boundSP.name))
                .collect(Collectors.toList());

        // for each unbound search problem filter if all positions have only one amino acid
        List<SearchProblem> mutableUnboundSp = unBoundSearchProblems.stream()
                .filter(sp -> getAllowedAA(sp).stream().allMatch(pos -> pos.size() == 1))
                .collect(Collectors.toList());

        // make sure we only have one mutable strand
        if (mutableUnboundSp.size() != 1) {
            throw new RuntimeException("ERROR: Only One Strand Can Currently Be NonMutable");
        }
        return mutableUnboundSp.get(0);
    }

    private ArrayList<ArrayList<String>> getAllowedAA(SearchProblem sp) {
        //First get all allowed AA's for bound position
        ArrayList<ArrayList<String>> allowedAAsComplex = cfp.getAllowedAAs();
        ArrayList<String> flexResNumComplex = cfp.getFlexRes();

        SearchProblem boundSP = getBoundSearchProblem(this.searchSpaces);

        //Add WT AA's if we need to
        for (int posNum = 0; posNum < boundSP.confSpace.numPos; posNum++) {
            ArrayList<String> currentAAOptions = allowedAAsComplex.get(posNum);
            if (cfp.params.getBool("AddWT")) {
                Residue res = boundSP.confSpace.m.getResByPDBResNumber(flexResNumComplex.get(posNum));
                if (!StringParsing.containsIgnoreCase(allowedAAsComplex.get(posNum), res.template.name)) {
                    currentAAOptions.add(res.template.name);
                }
            }
            if (currentAAOptions.isEmpty()) {
                throw new RuntimeException("ERROR: No AAtype for Position: " + posNum);
            }
        }
        int numPosBound = allowedAAsComplex.size();

        // Now depending on the search problem specified, get the appropriate allowed AAs
        ArrayList<ArrayList<String>> allowedAAs;

        boolean isBound = sp.name.equals(boundSP.name);
        //is the search problem specified the bound search problem?
        if (isBound) {
            allowedAAs = allowedAAsComplex;
        } else {
            // if not bound, is it strand0 or strand1?
            // unbound strand names should end with strand0 or strand1
            boolean isStrand0 = sp.name.endsWith("strand0");
            if (!(isStrand0 || sp.name.endsWith("strand1"))) {
                throw new RuntimeException("Error: Strand Naming Conventions Appear to Be Off");
            }

            int numFlex = sp.confSpace.numPos;
            if (isStrand0) {
                allowedAAs = new ArrayList(allowedAAsComplex.subList(0, numFlex));
            } else {
                allowedAAs = new ArrayList(allowedAAsComplex.subList(numPosBound - numFlex, numPosBound));
            }
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

    //Make sure each mutable state has the same Allowed AA
    //Return Allowed AA for all mutable positions
    private ArrayList<ArrayList<String>> handleAATypeOptions(ArrayList<ArrayList<String>> boundAllowedAAs, ArrayList<ArrayList<String>> unboundAllowedAAs) {
        ArrayList<ArrayList<String>> AATypeOptionsBound = boundAllowedAAs.stream().filter(aaTypes -> aaTypes.size() > 1).collect(Collectors.toCollection(ArrayList::new));

        ArrayList<ArrayList<String>> AATypeOptionsUnbound = unboundAllowedAAs.stream().filter(aaTypes -> aaTypes.size() > 1).collect(Collectors.toCollection(ArrayList::new));
        if (AATypeOptionsUnbound.size() != AATypeOptionsBound.size()) {
            throw new RuntimeException("ERROR: Different Number of Mutable Positions between Bound and Unbound");
        }
        for (int posNum = 0; posNum < AATypeOptionsUnbound.size(); posNum++) {
            ArrayList<String> AATypesForPos = AATypeOptionsUnbound.get(posNum);
            for (int aaIndex = 0; aaIndex < AATypesForPos.size(); aaIndex++) {
                if (!(AATypeOptionsBound.get(posNum).contains(AATypesForPos.get(aaIndex)))) {
                    throw new RuntimeException("ERROR: AAType Different for Bound and Unbound Mutable Residues");
                }
            }
        }
        
        return AATypeOptionsBound;
    }

    //Loads energy matrices and prune for COMETS
    private void loadEMatandPrune(double pruningInterval) {
        cfp.params.setValue("TYPEDEP", "TRUE");
        for (int state = 0; state < searchSpaces.length; state++) {
            SearchProblem searchProblem = this.searchSpaces[state];

            System.out.println("Precomputing Energy Matrix for " + searchProblem.name + " state");
            searchProblem.loadEnergyMatrix();
            System.out.println("Initializing Pruning for " + searchProblem.name + " state");
            initializePruning(searchProblem);
            PruningControl pruning = cfp.setupPruning(searchProblem, pruningInterval, false, false);
            pruning.prune();
        }
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

    //For testing/debugging purposes
    public void exhaustiveSublinearKStarSearch() {

        System.out.println();
        System.out.println("CHECKING Sulinear KStar RESULT BY EXHAUSTIVE SEARCH");
        System.out.println();

        ArrayList<String[]> seqList = listAllSeqs();
        int numSeqs = seqList.size();
        double stateScore[][] = new double[numSeqs][2];

        for (int seqNum = 0; seqNum < numSeqs; seqNum++) {
            String[] sequence = seqList.get(seqNum);
            for (String aa : sequence) {
                System.out.print(aa + " ");
            }
            for (int state = 0; state < this.numStates; state++) {
                SearchProblem searchProb = state == 0 ? boundSearchSpace : unboundMutSearchSpace;

                ArrayList<Integer> posNumsForState = new ArrayList<>();
                for (int pos = 0; pos < searchProb.confSpace.numPos; pos++) {
                    posNumsForState.add(pos);
                }

                PruningMatrix seqPruneMat = new PruningMatrix(searchProb.pruneMat.getSubsetMatrix(posNumsForState));
                handleStatePruning(seqPruneMat, searchProb.confSpace, sequence, this.mutable2StatePosNums.get(state));

                boolean allRotsPruned = false;
                for (int pos = 0; pos < seqPruneMat.numPos(); pos++) {
                    if (seqPruneMat.unprunedRCsAtPos(pos).isEmpty()) {
                        allRotsPruned = true;
                        System.out.println("All Rots Pruned");
                    }
                }
                if (allRotsPruned) {
                    stateScore[seqNum][state] = Double.POSITIVE_INFINITY;
                } else {
                    TRBP.setNumEdgeProbUpdates(3);
                    PartFuncTree tree = new PartFuncTree(searchProb.emat, seqPruneMat);
                    stateScore[seqNum][state] = -tree.computeEpsilonApprox(0.5);
                }
            }
            System.out.print("  Score: " + (stateScore[seqNum][0] - stateScore[seqNum][1] + objFcn.constTerm) + "\n");
        }

        //now find the best sequence and obj fcn value
        int topSeqNum = -1;
        double bestSeqScore = Double.POSITIVE_INFINITY;

        for (int seqNum = 0;
                seqNum < numSeqs;
                seqNum++) {
            boolean constrSatisfied = true;
            for (LME constr : gmecConstraints) {
                if (constr.eval(stateScore[seqNum]) > 0) {
                    constrSatisfied = false;
                }
            }

            if (constrSatisfied) {
                double seqScore = objFcn.eval(stateScore[seqNum]);
                if (seqScore < bestSeqScore) {
                    bestSeqScore = seqScore;
                    topSeqNum = seqNum;
                }
            }
        }
        this.bestSequence = seqList.get(topSeqNum);
        System.out.println();
        if (topSeqNum
                == -1) {
            System.out.println("EXHAUSTIVE SEARCH FINDS NO CONSTR-SATISFYING SEQUENCES");
        } else {
            System.out.println("EXHAUSTIVE SublinearKStar BEST SCORE: " + bestSeqScore + " SEQUENCE: ");
            for (String aa : seqList.get(topSeqNum)) {
                System.out.println(aa);
            }
        }

        System.out.println();
    }

    private void handleStatePruning(PruningMatrix pruneMat, ConfSpace cSpace, String[] sequence, ArrayList<Integer> mutable2StatePosNum) {
        if (mutable2StatePosNum.size() != sequence.length) {
            throw new RuntimeException("ERROR: Length of Sequence Not Equal to NumMutablePositions. Cannot complete exhaustive search");
        }

        for (int i = 0; i < mutable2StatePosNum.size(); i++) {
            int posNum = mutable2StatePosNum.get(i);
            String AAtype = sequence[i];
            for (int rc : pruneMat.unprunedRCsAtPos(posNum)) {
                String rcAAType = cSpace.posFlex.get(posNum).RCs.get(rc).AAType;

                if (!rcAAType.equalsIgnoreCase(AAtype)) {
                    pruneMat.markAsPruned(new RCTuple(posNum, rc));
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
}
