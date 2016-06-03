/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.AStarTree;
import edu.duke.cs.osprey.astar.comets.LME;
import edu.duke.cs.osprey.astar.comets.UpdatedPruningMatrix;
import edu.duke.cs.osprey.astar.partfunc.PartFuncTree;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.pruning.Pruner;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;
import org.apache.commons.lang.ArrayUtils;

/**
 *
 * @author hmn5
 */
public class SublinearKStarTree extends AStarTree {

    public int numTreeLevels; //number of mutable residues

    LME objFcn; //objective function to minimize
    LME[] constraints; //constraints on our sequence

    public ArrayList<ArrayList<String>> AATypeOptions; //The allowed amino-acids at each level

    int numMaxMut; //number of mutatations allowed away from wtSeq (-1 means no cap)
    String wtSeq[]; //wt sequence
    //information on states
    int numStates;//how many states there are
    //they have to have the same mutable residues & AA options,
    //though the residues involved may be otherwise different

    public SearchProblem[] mutableSearchProblems;//SearchProblems involved in COMETS search

    public SearchProblem nonMutableSearchProblem;

    ArrayList<ArrayList<Integer>> mutable2StatePosNums;
    //mutable2StatePosNum.get(state) maps levels in this tree to flexible positions for state
    //(not necessarily an onto mapping)

    int stateNumPos[];

    //Maps the bound res num to the corresponding unbound emat
    HashMap<Integer, EnergyMatrix> boundResNumToUnboundEmat;
    //Maps the bound res num to the corresponding unbound res num
    HashMap<Integer, Integer> boundPosNumToUnboundPosNum;
    //Maps the bound res num to boolean that is true if res num is part of mutable
    //strand
    HashMap<Integer, Boolean> boundResNumToIsMutableStrand;
    //determines if two residues are on the same strand
    boolean[][] belongToSameStrand;

    boolean useComets = false;

    int numSamplesGumbel = 50;
    final double constRT = PoissonBoltzmannEnergy.constRT;

    double bestChangeFreeEnergy;
    String[] bestSequence;

    double pruningInterval = 50;

    public int numLeafNodesVisited = 0;
    public int numLeafNodesExpanded = 0;

    public SublinearKStarTree(int numTreeLevels, LME objFcn, LME[] constraints,
            ArrayList<ArrayList<String>> AATypeOptions, int numMaxMut, String[] wtSeq,
            int numStates, SearchProblem[] stateSP, SearchProblem nonMutableSearchProblem,
            ArrayList<ArrayList<Integer>> mutable2StatePosNums) {

        this.numTreeLevels = numTreeLevels;
        this.objFcn = objFcn;
        this.constraints = constraints;
        this.AATypeOptions = AATypeOptions;
        this.numMaxMut = numMaxMut;
        this.wtSeq = wtSeq;
        this.numStates = numStates;
        this.mutableSearchProblems = stateSP;
        this.nonMutableSearchProblem = nonMutableSearchProblem;
        this.mutable2StatePosNums = mutable2StatePosNums;

        stateNumPos = new int[numStates];
        for (int state = 0; state < numStates; state++) {
            stateNumPos[state] = stateSP[state].confSpace.numPos;
        }

        this.boundResNumToUnboundEmat = getBoundPosNumToUnboundEmat();
        this.boundPosNumToUnboundPosNum = getBoundPosNumToUnboundPosNum();
        this.boundResNumToIsMutableStrand = getBoundPosNumberToIsMutableStrand();
        this.belongToSameStrand = getSameStrandMatrix();

    }

    private double boundFreeEnergyChange(SequenceNode seqNode) {
        if (seqNode.isFullyDefined())//fully-defined sequence
        {
            double score = 0;
            for (int state = 0; state < this.numStates; state++) {
                if (seqNode.stateTrees[state] == null) {
                    return Double.POSITIVE_INFINITY;
                }
            }
            score += seqNode.stateLBFreeEnergy[0];
            score -= seqNode.stateUBFreeEnergy[1];
            if (score == 0) {
                System.out.println("Score = 0");
            }
            return score;
        } else {
//            double realBound = computeExactBoundPerSequence(seqNode);
//            double maxIntBound = computeMaxInterfacePerSequence(seqNode, mutableSearchProblems[0].emat, seqNode.pruneMat[0]);
            double logKStarUB = calcMaxInterfaceScore(seqNode);
//            return logKStarUB;
            return logKStarUB;
        }
    }

    @Override
    public ArrayList<AStarNode> getChildren(AStarNode curNode) {
        SequenceNode seqNode = (SequenceNode) curNode;
        ArrayList<AStarNode> ans = new ArrayList<>();

        if (seqNode.isFullyDefined()) {
            seqNode.expandConfTree();
            seqNode.setScore(boundFreeEnergyChange(seqNode));
            ans.add(seqNode);
            this.numLeafNodesExpanded++;
            return ans;
        } else {
            //expand next position...
            int[] curAssignments = seqNode.getNodeAssignments();

            for (int splitPos = 0; splitPos < numTreeLevels; splitPos++) {
                if (curAssignments[splitPos] < 0) {//we can split this level

                    for (int aa = 0; aa < AATypeOptions.get(splitPos).size(); aa++) {
                        int[] childAssignments = curAssignments.clone();
                        childAssignments[splitPos] = aa;
                        UpdatedPruningMatrix[] childPruneMat = new UpdatedPruningMatrix[numStates];
                        for (int state = 0; state < numStates; state++) {
                            childPruneMat[state] = doChildPruning(state, seqNode.pruneMats[state], splitPos, aa);
                        }

                        SequenceNode childNode = new SequenceNode(childAssignments, childPruneMat);

                        if (splitPos == numTreeLevels - 1) {//sequence is fully defined...make conf trees
                            makeSeqPartFuncTree(childNode);
                             this.numLeafNodesVisited++;
                        }
                        printSequence(getSequence(childNode));
                        childNode.setScore(boundFreeEnergyChange(childNode));
                        System.out.println("Score: " + childNode.getScore());
                        ans.add(childNode);
                    }

                    return ans;
                }
            }
            throw new RuntimeException("ERROR: Not splittable position found but sequence not fully defined...");
        }
    }

    private double calcSequenceScore(SequenceNode seqNode) {
        if (seqNode.stateTrees[0] == null || seqNode.stateTrees[1] == null) {
            return Double.POSITIVE_INFINITY;
        }

        SearchProblem boundSP = this.mutableSearchProblems[0];
        SearchProblem ligandSP = this.mutableSearchProblems[1];

        EnergyMatrix boundEmat = boundSP.useTupExpForSearch ? boundSP.tupExpEMat : boundSP.emat;
        EnergyMatrix unboundEmat = ligandSP.useTupExpForSearch ? ligandSP.tupExpEMat : ligandSP.emat;

//            double Ebound = computeLogZGumbel(boundEmat, seqNode.pruneMat[0]);
//            double EunBound = computeLogZGumbel(unboundEmat, seqNode.pruneMat[1]);
        PartFuncTree treeB = new PartFuncTree(boundEmat, seqNode.pruneMats[0]);
        PartFuncTree treeU = new PartFuncTree(unboundEmat, seqNode.pruneMats[1]);
        double Ebound = treeB.computeEpsilonApprox(0.1);
        double EunBound = treeU.computeEpsilonApprox(0.1);
        double score = Ebound - EunBound + objFcn.getConstTerm();

        return -score;

    }

    private double calcMaxInterfaceScore(SequenceNode seqNode) {
        SearchProblem boundSP = mutableSearchProblems[0];
        SearchProblem ligandSP = mutableSearchProblems[1];

        ArrayList<Integer> boundPosNums = getAllBoundPosNums();
        ArrayList<Integer> proteinBoundPosNums = getProteinPosNums(true);
        ArrayList<Integer> ligandBoundPosNums = getLigandPosNums(true);

        boolean[][] interactionGraph_protein = createInteractionGraph(boundPosNums, proteinBoundPosNums, proteinBoundPosNums);
        boolean[][] interactionGraph_p_l = createInteractionGraph(boundPosNums, proteinBoundPosNums, ligandBoundPosNums);
        boolean[][] interactionGraph = addInteractionGraphs(interactionGraph_protein, interactionGraph_p_l);

        EnergyMatrix ematSubset;
        if (boundSP.useTupExpForSearch) {
            ematSubset = new EnergyMatrix(boundSP.tupExpEMat.getSubsetMatrix(boundPosNums));
            ematSubset.updateMatrixCrossTerms(interactionGraph);
            ematSubset.addInternalEnergies(boundSP.tupExpEMat, proteinBoundPosNums);
            ematSubset.addCrossTermInternalEnergies(boundSP.tupExpEMat, ligandSP.tupExpEMat, ligandBoundPosNums, boundPosNumToUnboundPosNum);
        } else {
            ematSubset = new EnergyMatrix(boundSP.emat.getSubsetMatrix(boundPosNums));
            ematSubset.updateMatrixCrossTerms(interactionGraph);
            ematSubset.addInternalEnergies(boundSP.emat, proteinBoundPosNums);
            ematSubset.addCrossTermInternalEnergies(boundSP.emat, ligandSP.emat, ligandBoundPosNums, boundPosNumToUnboundPosNum);
        }
//        PruningMatrix pruneMatSubset = new PruningMatrix(seqNode.pruneMats[0].getSubsetMatrix(boundPosNums));

        double score;
//            PartFuncTree tree = new PartFuncTree(ematSubset, pruneMatSubset);
//            score = -computeLogZGumbel(ematSubset, pruneMatSubset);
//            score = -tree.computeEpsilonApprox(0.1);
        TRBP.setNumEdgeProbUpdates(0);
//        System.out.println("Unpruned: " + seqNode.pruneMats[0].unprunedRCsAtPos(5));


        /*        ConfTree tree = new ConfTree(ematSubset, seqNode.pruneMats[0]);
         int[] gmec = tree.nextConf();
         if (gmec != null) {
         RCTuple tup = new RCTuple(gmec);
         System.out.println("GMEC: " + ematSubset.getInternalEnergy(tup) / this.constRT);
         String aa = boundSP.confSpace.posFlex.get(7).RCs.get(gmec[7]).AAType;
         boundSP.outputMinimizedStruct(gmec, aa + ".x.pdb");
         } else {
         System.out.println("GMEC NULL");
         }
         */
//        ConfTree tree2 = new ConfTree(boundSP.emat, seqNode.pruneMats[0]);
//        int[] gmec2 = tree2.nextConf();
//        RCTuple tup2 = new RCTuple(gmec2);
//.        System.out.println("GMEC: " + boundSP.emat.getInternalEnergy(tup2) / this.constRT);
//        boundSP.outputMinimizedStruct(gmec2, "trp.gmec.pdb");
//        if (boundSP.emat.getInternalEnergy(tup2) > 100) {
//            throw new RuntimeException();
//        }
        if (hasUnprunedRCs(seqNode.pruneMats[0])) {
            MarkovRandomField mrf = new MarkovRandomField(ematSubset, seqNode.pruneMats[0], 0.0);
            TRBP trbp = new TRBP(mrf);
            score = -trbp.getLogZ();
        } else {
            score = Double.POSITIVE_INFINITY;
        }
        return score + objFcn.getConstTerm();
    }

    private boolean hasUnprunedRCs(PruningMatrix pm) {
        for (int pos = 0; pos < pm.numPos(); pos++) {
            if (pm.unprunedRCsAtPos(pos).isEmpty()) {
                return false;
            }
        }
        return true;
    }

    private double computeExactBoundPerSequence(SequenceNode seqNode) {
        int[][] seqList = getAllSequences(seqNode);
        double score = Double.POSITIVE_INFINITY;
        for (int[] seq : seqList) {
            UpdatedPruningMatrix pruneMatSeqBound = new UpdatedPruningMatrix(seqNode.pruneMats[0]);
            UpdatedPruningMatrix pruneMatSeqUnbound = new UpdatedPruningMatrix(seqNode.pruneMats[1]);
            updateBoundPruneMatAtSequence(seqNode, pruneMatSeqBound, seq);
            updateUnoundPruneMatAtSequence(seqNode, pruneMatSeqUnbound, seq);
            PartFuncTree pfTreeBound = new PartFuncTree(mutableSearchProblems[0].emat, pruneMatSeqBound);
            PartFuncTree pfTreeUnbound = new PartFuncTree(mutableSearchProblems[1].emat, pruneMatSeqUnbound);

            double boundScore = -pfTreeBound.computeEpsilonApprox(0.1);
            double unboundScore = -pfTreeUnbound.computeEpsilonApprox(0.1);
//            double boundScore = -computeLogZGumbel(mutableSearchProblems[0].emat, pruneMatSeqBound);
//            double unboundScore = -computeLogZGumbel(mutableSearchProblems[1].emat, pruneMatSeqUnbound);
            System.out.print("Sequence: ");
            for (int pos = 0; pos < seq.length; pos++) {
                System.out.print(this.AATypeOptions.get(pos).get(seq[pos]) + " ");
            }
            System.out.println("   Score: " + (boundScore - unboundScore));
//            DiscretePartFunc dpf = new DiscretePartFunc(mutableSearchProblems[0].emat, pruneMatSeqBound, 0.1);
//            double exactScoreBound = -dpf.getLogZ();
//            DiscretePartFuncCalc dpfBound = new DiscretePartFuncCalc(mutableSearchProblems[0].emat, pruneMatSeqBound, 100, 0.1);
//            double exactScoreBound = dpfBound.calcLogZ();
//            DiscretePartFuncCalc dpfUnbound = new DiscretePartFuncCalc(mutableSearchProblems[1].emat, pruneMatSeqUnbound, 100, 0.1);
//            double exactScoreUnbound = dpfUnbound.calcLogZ();
//            DiscretePartFunc dfpU = new DiscretePartFunc(mutableSearchProblems[1].emat, pruneMatSeqUnbound, 0.1);
//            double exactScoreUnbound = -dfpU.getLogZ();
//            System.out.println("Exact Score: "+(exactScoreBound-exactScoreUnbound));
            score = Math.min(score, (boundScore - unboundScore));
        }
        return score;
    }

    /**
     * Creates an interaction graph over the set of residues defined by
     * allPositions Every element in subsetI is interacting with every element
     * in subsetJ interactions = {(i,j) | i in subsetI and j in subsetJ}
     *
     * @param allPositions all the positions that we are considering sorted by
     * original posNumber
     * @param subsetI list of positions interacting with subsetJ
     * @param subsetJ list of positions interacting with subsetI
     * @return interaction graph interactions = {(i,j) | i in subsetI and j in
     * subsetJ}
     */
    private boolean[][] createInteractionGraph(ArrayList<Integer> allPositions, ArrayList<Integer> subsetI, ArrayList<Integer> subsetJ) {
        int numPos = allPositions.size();
        boolean[][] interactionGraph = new boolean[numPos][numPos];

        //Initialize interactoin graph
        for (int posI = 0; posI < numPos; posI++) {
            for (int posJ = 0; posJ < numPos; posJ++) {
                interactionGraph[posI][posJ] = false;
            }
        }

        for (int posI : subsetI) {
            int newPosNum_I = allPositions.indexOf(posI);
            for (int posJ : subsetJ) {
                //We don't consider self-interactions here so we will leave it as false
                if (posI != posJ) {
                    int newPosNum_J = allPositions.indexOf(posJ);
                    interactionGraph[newPosNum_I][newPosNum_J] = true;
                    interactionGraph[newPosNum_J][newPosNum_I] = true;
                }
            }
        }
        return interactionGraph;
    }

    /**
     * Given two interaction graphs over the same graph we "add" them By this I
     * mean, (i,j) is interacting if (i,j) is interacting in either graphI or
     * graphJ
     *
     * @param interactionGraphI
     * @param interactionGraphJ
     * @return
     */
    private boolean[][] addInteractionGraphs(boolean[][] interactionGraphI, boolean[][] interactionGraphJ) {
        //Check to make sure each graph is the same size
        if (interactionGraphI.length != interactionGraphJ.length) {
            throw new RuntimeException("ERROR: Cannot add two interaction graphs of different size");
        }

        int numPos = interactionGraphI.length;
        boolean[][] interactionGraph = new boolean[numPos][numPos];
        for (int i = 0; i < numPos; i++) {
            for (int j = i + 1; j < numPos; j++) {
                interactionGraph[i][j] = (interactionGraphI[i][j] || interactionGraphJ[i][j]);
                interactionGraph[j][i] = interactionGraph[i][j];
            }
        }
        return interactionGraph;
    }

    private void makeSeqPartFuncTree(SequenceNode node) {
        //Given a node with a fully defined sequence, build its conformational search trees
        //for each state
        //If a state has no viable conformations, leave it null, with stateUB[state] = inf

        node.stateTrees = new PartFuncTree[numStates];
        node.stateUBFreeEnergy = new double[numStates];
        node.stateLBFreeEnergy = new double[numStates];

        for (int state = 0; state < numStates; state++) {

            //first make sure there are RCs available at each position
            boolean RCsAvailable = true;
            for (int pos = 0; pos < stateNumPos[state]; pos++) {
                if (node.pruneMats[state].unprunedRCsAtPos(pos).isEmpty()) {
                    RCsAvailable = false;
                    break;
                }
            }

            if (RCsAvailable) {
                TRBP.setNumEdgeProbUpdates(3);
                node.stateTrees[state] = new PartFuncTree(mutableSearchProblems[state].emat, node.pruneMats[state]);

                AStarNode rootNode = node.stateTrees[state].rootNode();
                node.stateTrees[state].initQueue(rootNode);//allocate queue and add root node
                node.stateLBFreeEnergy[state] = -node.stateTrees[state].getCurrentUpperBoundLogZ();
                node.stateUBFreeEnergy[state] = -node.stateTrees[state].getCurrentLowerBoundLogZ();
            } else {//no confs available for this state!
                node.stateTrees[state] = null;
                node.stateUBFreeEnergy[state] = Double.POSITIVE_INFINITY;
                node.stateLBFreeEnergy[state] = Double.NEGATIVE_INFINITY;
            }
        }
    }

    /**
     * Prunes rotamers to reflex the new allowed amino-acids
     *
     * @param state bound vs unbound state
     * @param parentMat parent pruning matrix
     * @param splitPos position that was split
     * @param aa amino acid label for new positions that was split
     * @return updated pruning matrix
     */
    private UpdatedPruningMatrix doChildPruning(int state, PruningMatrix parentMat, int splitPos, int aa) {
        //Create an update to parentMat (without changing parentMat)
        //to reflect that splitPos has been assigned an amino-acid type

        String assignedAAType = AATypeOptions.get(splitPos).get(aa);

        UpdatedPruningMatrix ans = new UpdatedPruningMatrix(parentMat);
        int posAtState = mutable2StatePosNums.get(state).get(splitPos);

        //first, prune all other AA types at splitPos
        for (int rc : parentMat.unprunedRCsAtPos(posAtState)) {
            //HUNTER: TODO: AATYperPerRes should only be one residue for now
            //We should change this to allow for rcs at the sequence search level
            String rcAAType = mutableSearchProblems[state].confSpace.posFlex.get(posAtState).RCs.get(rc).AAType;

            if (!rcAAType.equalsIgnoreCase(assignedAAType)) {
                ans.markAsPruned(new RCTuple(posAtState, rc));
            }
        }

        //now do any consequent singles & pairs pruning
        int numUpdates = ans.countUpdates();
        int oldNumUpdates;

        Pruner dee = new Pruner(mutableSearchProblems[state], ans, true, Double.POSITIVE_INFINITY,
                this.pruningInterval, false, mutableSearchProblems[state].useTupExpForSearch, false);
        //this is rigid, type-dependent pruning aiming for sequence GMECs
        //So Ew = Ival = 0

        do {//repeat as long as we're pruning things
            oldNumUpdates = numUpdates;
            dee.prune("GOLDSTEIN");
//            dee.prune("GOLDSTEIN PAIRS FULL");
            numUpdates = ans.countUpdates();
        } while (numUpdates > oldNumUpdates);

        return ans;

    }

    @Override
    public boolean isFullyAssigned(AStarNode node) {
        //HMN: TODO: We need to decide when a K* score is fully calculated

        if (!node.isFullyDefined()) {
            return false;
        }

        double lowerBound = this.pq.peek().getScore();
        SequenceNode seqNode = (SequenceNode) node;
        if (seqNode.stateUBFreeEnergy[0] - seqNode.stateLBFreeEnergy[1] < lowerBound) {
            this.bestChangeFreeEnergy = seqNode.getScore();
            this.bestSequence = getSequence(node);

            System.out.println("GMEC FOUND");
            return true;
        }

        double ub = seqNode.stateUBFreeEnergy[0] - seqNode.stateLBFreeEnergy[1];
        System.out.println("Upper Bound: " + ub);
        System.out.println("Lower Bound: " + lowerBound);
        System.out.println();

        if (this.pq.isEmpty()) {
            throw new RuntimeException("Queue is empty");
        }
        //For Ensembles, a fully-defined node will have an exact K* score for now...
        return false;
    }

    String[] getSequence(AStarNode node) {
        int[] assignment = node.getNodeAssignments();
        String[] seq = new String[assignment.length];
        for (int i = 0; i < assignment.length; i++) {
            seq[i] = this.AATypeOptions.get(i).get(assignment[i]);
        }
        return seq;
    }

    @Override
    public AStarNode rootNode() {
        int[] conf = new int[numTreeLevels];
        Arrays.fill(conf, -1);//indicates the sequence is not assigned

        PruningMatrix[] pruneMats = new PruningMatrix[numStates];
        for (int state = 0; state < numStates; state++) {
            pruneMats[state] = mutableSearchProblems[state].pruneMat;
        }

        SequenceNode root = new SequenceNode(conf, pruneMats);
        return root;
    }

    @Override
    public boolean canPruneNode(AStarNode node) {
        //Check if node can be pruned
        //This is traditionally based on constraints, thought we could pruned nodes
        //that are provably suboptimal

        //TODO: Implement constraints as in COMETS if desired
        SequenceNode seqNode = (SequenceNode) node;

        if (numMaxMut != -1) {
            //cap on number of mutations
            int mutCount = 0;
            int assignments[] = seqNode.getNodeAssignments();

            for (int level = 0; level < numTreeLevels; level++) {
                if (assignments[level] >= 0) {//AA type at level is assigned
                    if (!AATypeOptions.get(level).get(assignments[level]).equalsIgnoreCase(wtSeq[level]))//and is different from wtSeq
                    {
                        mutCount++;
                    }
                }
            }

            if (mutCount > numMaxMut)//prune based on number of mutations
            {
                return true;
            }
        }
        if (Double.isInfinite(seqNode.getScore())) {
            return true;
        }
        //TODO: for (LME constr : constraints) {....

        return false;
    }

    /**
     * For flexible position in bound matrix (0,1,2,...,numRes-1) we map to the
     * corresponding unbound energy matrix
     *
     * @return
     */
    public HashMap<Integer, EnergyMatrix> getBoundPosNumToUnboundEmat() {
        SearchProblem boundState = this.mutableSearchProblems[0];
        //Get res number from each flexible position in the bound state
        List<String> resNumsBound = boundState.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound mutable state
        SearchProblem unBoundMutableState = this.mutableSearchProblems[1];
        List<String> resNumsUnboundMutable = unBoundMutableState.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound non-mutable state
        SearchProblem unBoundNonMutableState = this.nonMutableSearchProblem;
        List<String> resNumsUnboundNonMutable = unBoundNonMutableState.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));
        //Map to corresponding energy matrix
        List<EnergyMatrix> unboundEmatPerPos = resNumsBound.stream()
                .map(posNum -> ArrayUtils.contains(resNumsUnboundMutable.toArray(), posNum) ? unBoundMutableState.emat : unBoundNonMutableState.emat)
                .collect(Collectors.toCollection(ArrayList::new));
        HashMap<Integer, EnergyMatrix> boundPosNumToUnboundEmat = new HashMap<>();
        for (int posNum = 0; posNum < unboundEmatPerPos.size(); posNum++) {
            boundPosNumToUnboundEmat.put(posNum, unboundEmatPerPos.get(posNum));
        }
        return boundPosNumToUnboundEmat;
    }

    /**
     * For flexible position in bound matrix (0,1,2,...,numRes-1) we map to the
     * corresponding position number in the unbound matrix
     *
     * @return
     */
    public HashMap<Integer, Integer> getBoundPosNumToUnboundPosNum() {
        SearchProblem boundState = this.mutableSearchProblems[0];
        //Get res number from each flexible position in the bound state
        List<String> resNumsBound = boundState.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound mutable state
        SearchProblem unBoundMutableState = this.mutableSearchProblems[1];
        List<String> resNumsUnboundMutable = unBoundMutableState.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound non-mutable state
        SearchProblem unBoundNonMutableState = this.nonMutableSearchProblem;
        List<String> resNumsUnboundNonMutable = unBoundNonMutableState.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));
        //Map to corresponding unbound position number
        List<Integer> unboundPosNumsPerPos = resNumsBound.stream()
                .map(posNum -> ArrayUtils.contains(resNumsUnboundMutable.toArray(), posNum) ? ArrayUtils.indexOf(resNumsUnboundMutable.toArray(), posNum)
                                : ArrayUtils.indexOf(resNumsUnboundNonMutable.toArray(), posNum))
                .collect(Collectors.toCollection(ArrayList::new));

        HashMap<Integer, Integer> boundPosNumToUnboundPosNum = new HashMap<>();
        for (int posNum = 0; posNum < unboundPosNumsPerPos.size(); posNum++) {
            boundPosNumToUnboundPosNum.put(posNum, unboundPosNumsPerPos.get(posNum));
        }
        return boundPosNumToUnboundPosNum;
    }

    /**
     * For flexible position in bound matrix (0,1,2,...,numRes-1) we map to the
     * corresponding strand number
     *
     * @return
     */
    public boolean[][] getSameStrandMatrix() {
        SearchProblem boundState = this.mutableSearchProblems[0];
        //Get res number from each flexible position in the bound state
        List<String> resNumsBound = boundState.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toList());

        //Get res number for each flexible position in the unbound mutable state
        SearchProblem unBoundMutableState = this.mutableSearchProblems[1];
        List<String> resNumsUnboundMutable = unBoundMutableState.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toList());

        //Get res number for each flexible position in the unbound non-mutable state
        SearchProblem unBoundNonMutableState = this.nonMutableSearchProblem;
        List<String> resNumsUnboundNonMutable = unBoundNonMutableState.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toList());
        //Map to corresponding unbound position number
        List<Integer> boundPosNumToStrandNum = resNumsBound.stream()
                .map(posNum -> ArrayUtils.contains(resNumsUnboundMutable.toArray(), posNum) ? 0 : 1)
                .collect(Collectors.toList());

        int numResidues = boundPosNumToStrandNum.size();
        boolean[][] belongToSameStrand = new boolean[numResidues][numResidues];
        for (int i = 0; i < numResidues; i++) {
            for (int j = i + 1; j < numResidues; j++) {
                if (boundPosNumToStrandNum.get(i).equals(boundPosNumToStrandNum.get(j))) {
                    belongToSameStrand[i][j] = true;
                    belongToSameStrand[j][i] = true;
                } else {
                    belongToSameStrand[i][j] = false;
                    belongToSameStrand[j][i] = false;
                }
            }
        }
        return belongToSameStrand;
    }

    /**
     * For flexible position in bound matrix (0,1,2,...,numRes-1) we map to the
     * boolean that is true if the res is part of boolean strand and false
     * otherwise
     *
     * @return
     */
    public HashMap<Integer, Boolean> getBoundPosNumberToIsMutableStrand() {
        SearchProblem boundState = this.mutableSearchProblems[0];
        //Get res number from each flexible position in the bound state
        List<String> resNumsBound = boundState.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound mutable state
        SearchProblem unBoundMutableState = this.mutableSearchProblems[1];
        List<String> resNumsUnboundMutable = unBoundMutableState.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound non-mutable state
        SearchProblem unBoundNonMutableState = this.nonMutableSearchProblem;
        List<String> resNumsUnboundNonMutable = unBoundNonMutableState.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));
        //Map to corresponding unbound position number
        List<Boolean> boundPosNumIsMutable = resNumsBound.stream()
                .map(posNum -> ArrayUtils.contains(resNumsUnboundMutable.toArray(), posNum))
                .collect(Collectors.toCollection(ArrayList::new));

        HashMap<Integer, Boolean> boundPosNumToIsMutableStrand = new HashMap<>();
        for (int posNum = 0; posNum < boundPosNumIsMutable.size(); posNum++) {
            boundPosNumToIsMutableStrand.put(posNum, boundPosNumIsMutable.get(posNum));
        }
        return boundPosNumToIsMutableStrand;
    }

    /**
     * This returns the list of integers corresponding to the assigned position
     * numbers of the ligand in the confSpace This is useful for computing GMECs
     * over partial spaces.
     *
     * @param seqNode the current sequence node
     * @param useBoundState do we want to position numbers corresponding to
     * bound state (true) or unbound state (false)
     * @return
     */
    private ArrayList<Integer> getLigandAssignedPosNums(SequenceNode seqNode, boolean useBoundState) {
        int[] partialSeq = seqNode.getNodeAssignments();
        //Get the mapping between mutable position in sequence node assignment and
        //the corresponding position number in the confSpace for bound or unbound ligand
        ArrayList<Integer> mutable2PosNum = useBoundState ? this.mutable2StatePosNums.get(0) : this.mutable2StatePosNums.get(1);
        //Get the corresponding searchProblem for bound or unbound ligand
        SearchProblem searchProblem = useBoundState ? this.mutableSearchProblems[0] : this.mutableSearchProblems[1];

        ArrayList<Integer> ligandAssignedPosNums = new ArrayList<>();

        //Iterate over flexible position numbers 
        for (int posNum = 0; posNum < searchProblem.confSpace.numPos; posNum++) {

            //If we are using the bound state, we must check if the position belongs to 
            //Ligand or protein
            //If we are not using the bound state then it must belong to ligand
            if (this.boundResNumToIsMutableStrand.get(posNum) || !(useBoundState)) {
                //position is on mutable strand, implying it belongs to ligand

                //Now check if the posNum is mutable
                if (mutable2PosNum.contains(posNum)) {
                    //It is mutable so get the index of mutable position wrt sequence node assignment
                    int index = mutable2PosNum.indexOf(posNum);
                    //Now we check if this mutable position is assigned
                    if (partialSeq[index] >= 0) {
                        //It is assigned so add it to our list
                        ligandAssignedPosNums.add(posNum);
                    }
                } else {//if the posNum is NOT mutable then it is assigned already assigned
                    //add since it is assigned and on ligand
                    ligandAssignedPosNums.add(posNum);
                }
            }
        }
        return ligandAssignedPosNums;
    }

    /**
     * This returns the list of integers corresponding to the assigned position
     * numbers of the ligand in the confSpace This is useful for computing GMECs
     * over partial spaces.
     *
     * @param seqNode current sequence node in tree
     * @param useBoundState do we want the posNums to correspond to the bound
     * state or unbound state
     * @return
     */
    private ArrayList<Integer> getLigandUnassignedPosNums(SequenceNode seqNode, boolean useBoundState) {
        int[] partialSeq = seqNode.getNodeAssignments();

        //Get the mapping between mutable position in sequence node assignment and
        //the corresponding position number in the confSpace for bound or unbound ligand
        ArrayList<Integer> mutable2PosNum = useBoundState ? this.mutable2StatePosNums.get(0) : this.mutable2StatePosNums.get(1);
        //Get the corresponding searchProblem for bound or unbound ligand
        SearchProblem searchProblem = useBoundState ? this.mutableSearchProblems[0] : this.mutableSearchProblems[1];

        ArrayList<Integer> ligandUnassignedPosNums = new ArrayList<>();

        //Iterate over flexible position numbers
        for (int posNum = 0; posNum < searchProblem.confSpace.numPos; posNum++) {

            //If we are using the bound state, we must check if the position belongs to 
            //Ligand or protein
            //If we are not using the bound state then it must belong to ligand
            if (this.boundResNumToIsMutableStrand.get(posNum) || !(useBoundState)) {
                //Position is on the mutable (ligand) strand

                //Now check if the posNum is mutable
                if (mutable2PosNum.contains(posNum)) {
                    //It is mutable so get the index of the mutable position wrt sequence node assignment
                    int index = mutable2PosNum.indexOf(posNum);
                    //Now we check if this mutable position is unassigned
                    if (partialSeq[index] < 0) {
                        //It is unassigned so add to our list
                        ligandUnassignedPosNums.add(posNum);
                    }
                }
            }
        }
        return ligandUnassignedPosNums;
    }

    /**
     * Returns the list of posNums corresponding to the protein residues for
     * either the bound search space or unbound search space
     *
     * @param useBoundState if true, use bound search space
     * @return
     */
    private ArrayList<Integer> getProteinPosNums(boolean useBoundState) {
        SearchProblem searchSpace = useBoundState ? mutableSearchProblems[0] : nonMutableSearchProblem;

        ArrayList<Integer> proteinPosNums = new ArrayList<>();

        for (int pos = 0; pos < searchSpace.confSpace.numPos; pos++) {
            //If we are using the bound state we must check if posNum
            //belongs to protein or ligand
            if (useBoundState) {
                //Check if pos belongs to nonmutable (protein) strand
                if (!this.boundResNumToIsMutableStrand.get(pos)) {
                    proteinPosNums.add(pos);
                }
            } else {//If we are using the unbound state, every posNum is part of protein
                proteinPosNums.add(pos);
            }
        }
        return proteinPosNums;
    }

    /**
     * Returns the list of posNums corresponding to the ligand residues for
     * either the bound search space or the unbound searchSpace
     *
     * @param useBoundState if true, use bound search space
     * @return
     */
    private ArrayList<Integer> getLigandPosNums(boolean useBoundState) {
        SearchProblem searchSpace = useBoundState ? mutableSearchProblems[0] : mutableSearchProblems[1];

        ArrayList<Integer> ligandPosNums = new ArrayList<>();

        for (int pos = 0; pos < searchSpace.confSpace.numPos; pos++) {
            if (useBoundState) {
                if (this.boundResNumToIsMutableStrand.get(pos)) {
                    ligandPosNums.add(pos);
                }
            } else {
                ligandPosNums.add(pos);
            }
        }
        return ligandPosNums;
    }

    /**
     * Returns a list of all posNums (simply 0,1,...,numPosNums-1) for the bound
     * state
     *
     * @return
     */
    private ArrayList<Integer> getAllBoundPosNums() {
        SearchProblem searchSpace = mutableSearchProblems[0];
        ArrayList<Integer> allPosNums = new ArrayList<>();
        for (int pos = 0; pos < searchSpace.confSpace.numPos; pos++) {
            allPosNums.add(pos);
        }
        return allPosNums;
    }

    boolean updateUB(int state, AStarNode expNode, SequenceNode seqNode) {
        //Get an upper-bound on the node by a little FASTER run, generating UBConf
        //store UBConf and UB in expNode
        //expNode is in seqNode.stateTrees[state]
        //we'll start with the starting conf (likely from a parent) if provided
        //return true unless no valid conf is possible...then false

        int assignments[] = expNode.getNodeAssignments();
        ArrayList<ArrayList<Integer>> allowedRCs = new ArrayList<>();

        for (int pos = 0; pos < stateNumPos[state]; pos++) {

            ArrayList<Integer> posOptions = new ArrayList<>();

            if (assignments[pos] == -1) {
                posOptions = seqNode.pruneMats[state].unprunedRCsAtPos(pos);
                if (posOptions.isEmpty())//no options at this position
                {
                    return false;
                }
            } else//assigned
            {
                posOptions.add(assignments[pos]);
            }

            allowedRCs.add(posOptions);
        }

        int startingConf[] = expNode.UBConf;
        int[] UBConf = startingConf.clone();
        //ok first get rid of anything in startingConf not in expNode's conf space,
        //replace with the lowest-intra+shell-E conf
        for (int level = 0; level < stateNumPos[state]; level++) {

            if (!allowedRCs.get(level).contains(startingConf[level])) {
                //if( ! levelOptions.get(level).get(expNode.conf[level]).contains( startingConf[level] ) ){

                double bestISE = Double.POSITIVE_INFINITY;
                int bestRC = allowedRCs.get(level).get(0);
                for (int rc : allowedRCs.get(level)) {
                    double ise = getEnergyMatrix(state).getOneBody(level, rc);
                    if (ise < bestISE) {
                        bestISE = ise;
                        bestRC = rc;
                    }
                }

                UBConf[level] = bestRC;
            }
        }

        double curE = getEnergyMatrix(state).confE(UBConf);
        boolean done = false;

        while (!done) {

            done = true;

            for (int level = 0; level < stateNumPos[state]; level++) {

                int testConf[] = UBConf.clone();

                for (int rc : allowedRCs.get(level)) {
                    testConf[level] = rc;

                    //if(!canPrune(testConf)){//pruned conf unlikely to be good UB
                    //would only prune if using pruned pair flags in A*
                    double testE = getEnergyMatrix(state).confE(testConf);
                    if (testE < curE) {
                        curE = testE;
                        UBConf[level] = rc;
                        done = false;
                    }
                    //}
                }
            }
        }

        expNode.UBConf = UBConf;
        expNode.UB = getEnergyMatrix(state).confE(UBConf);

        return true;
    }

    private EnergyMatrix getEnergyMatrix(int state) {
        if (mutableSearchProblems[state].useTupExpForSearch)//Using LUTE
        {
            return mutableSearchProblems[state].tupExpEMat;
        } else//discrete flexibility w/o LUTE
        {
            return mutableSearchProblems[state].emat;
        }
    }

    void printBestSeqInfo(SequenceNode seqNode) {
        //About to return the given fully assigned sequence from A*
        //provide information
        System.out.println("SeqTree: A* returning conformation; lower bound = " + seqNode.getScore() + " nodes expanded: " + numExpanded);

        System.out.print("Sequence: ");

        for (int level = 0; level < numTreeLevels; level++) {
            System.out.print(AATypeOptions.get(level).get(seqNode.getNodeAssignments()[level]) + " ");
        }
        System.out.println();

        System.out.println();
        System.out.println(numExpanded + " expanded; " + getQueue().size() + " nodes in tree, of which "
                + numPruned + " pruned.");
    }

    @Override
    public int[] outputNode(AStarNode node) {
        System.out.println();
        printBestSeqInfo((SequenceNode) node);
        return node.getNodeAssignments();
    }

    public void printSequence(String[] sequence) {
        StringBuffer buffer = new StringBuffer();
        for (String aaType : sequence) {
            buffer.append(" " + aaType);
        }
        System.out.println(buffer);
    }

    public String[] getSequence(SequenceNode node) {
        int[] assignments = node.getNodeAssignments();
        int numMotPos = assignments.length;
        String[] sequence = new String[numMotPos];

        for (int mutPos = 0; mutPos < numMotPos; mutPos++) {
            int aaTypeVal = assignments[mutPos];
            String aaType;
            if (aaTypeVal == -1) {
                aaType = "XXX";
            } else {
                aaType = this.AATypeOptions.get(mutPos).get(aaTypeVal);
            }
            sequence[mutPos] = aaType;
        }
        return sequence;
    }

    int[][] getAllSequences(SequenceNode seqNode) {
        int[] assignmnet = seqNode.getNodeAssignments();

        //Iterate over mutable positions and get the total number of sequences
        int numSequences = 1;
        //Keep track of the last assigned level
        int lastAssignedLevel = -1;
        for (int i = 0; i < this.numTreeLevels; i++) {
            if (assignmnet[i] == -1) {
                numSequences = numSequences * this.AATypeOptions.get(i).size();
            } else {
                lastAssignedLevel = i;
            }
        }

        int[] currentSeq = ArrayUtils.subarray(assignmnet, 0, lastAssignedLevel + 1);
        int[][] seqList = new int[1][currentSeq.length];
        seqList[0] = currentSeq;

        for (int mutPos = lastAssignedLevel + 1; mutPos < this.numTreeLevels; mutPos++) {
            seqList = getAllSequencesHelper(seqList, mutPos);
        }
        return seqList;
    }

    int[][] getAllSequencesHelper(int[][] currentSeqList, int level) {
        //the number of current sequences
        int numCurrentSeqs = currentSeqList.length;
        //the number of new sequences after we add all AA's from mutable pos
        int numNewSeqs = numCurrentSeqs * this.AATypeOptions.get(level).size();
        //the current length of our sequences
        int currentSeqLength = currentSeqList[0].length;
        //the new length after we add mutPos
        int newSeqLength = currentSeqLength + 1;

        //the new array of sequences
        int[][] ans = new int[numNewSeqs][currentSeqList.length + 1];

        int sequenceNum = 0;
        for (int[] sequence : currentSeqList) {
            for (int aa = 0; aa < this.AATypeOptions.get(level).size(); aa++) {
                int[] newSequence = Arrays.copyOf(sequence, newSeqLength);
                //add AA from mutPos
                newSequence[newSeqLength - 1] = aa;
                ans[sequenceNum] = newSequence;
                //iterate
                sequenceNum++;
            }
        }

        return ans;
    }

    void updateBoundPruneMatAtSequence(SequenceNode seqNode, PruningMatrix pruneMat, int[] seq) {
        int[] assignments = seqNode.getNodeAssignments();
        SearchProblem boundSP = mutableSearchProblems[0];

        //Update the bound prune mat for this sequence
        for (int level = 0; level < assignments.length; level++) {
            if (assignments[level] == -1) {
                int posNum = mutable2StatePosNums.get(0).get(level);
                int aaNum = seq[level];
                String AAType = this.AATypeOptions.get(level).get(aaNum);
                for (int rot : pruneMat.unprunedRCsAtPos(posNum)) {
                    if (!boundSP.confSpace.posFlex.get(posNum).RCs.get(rot).AAType.equalsIgnoreCase(AAType)) {
                        pruneMat.markAsPruned(new RCTuple(posNum, rot));
                    }
                }
            }
        }
    }

    void updateUnoundPruneMatAtSequence(SequenceNode seqNode, PruningMatrix pruneMat, int[] seq) {
        int[] assignments = seqNode.getNodeAssignments();
        SearchProblem ligandSP = mutableSearchProblems[1];

        //Update the bound prune mat for this sequence
        for (int level = 0; level < assignments.length; level++) {
            if (assignments[level] == -1) {
                int posNum = mutable2StatePosNums.get(1).get(level);
                int aaNum = seq[level];
                String AAType = this.AATypeOptions.get(level).get(aaNum);
                for (int rot : pruneMat.unprunedRCsAtPos(posNum)) {
                    if (!ligandSP.confSpace.posFlex.get(posNum).RCs.get(rot).AAType.equalsIgnoreCase(AAType)) {
                        pruneMat.markAsPruned(new RCTuple(posNum, rot));
                    }
                }
            }
        }
    }

    public double getLogSequenceSpaceSize() {
        double logNumSeq = 0;
        for (int pos = 0; pos < this.AATypeOptions.size(); pos++) {
            logNumSeq += Math.log(this.AATypeOptions.get(pos).size());
        }
        return logNumSeq;
    }

}
