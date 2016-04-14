/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.comets;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.AStarTree;
import edu.duke.cs.osprey.astar.ConfTreeSuper;
import edu.duke.cs.osprey.astar.Mplp;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.SearchProblemSuper;
import edu.duke.cs.osprey.confspace.SuperRCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PrunerSuper;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.PriorityQueue;
import java.util.stream.Collectors;
import org.apache.commons.lang.ArrayUtils;

/**
 * Replica of COMETSTree but supports superRCs
 *
 * @author hmn5
 */
public class COMETSTreeSuper extends AStarTree {

    int numTreeLevels;//number of residues with sequence changes

    LME objFcn;//we are minimizing objFcn...
    LME[] constraints;//with the constraints that constr all <= 0

    ArrayList<ArrayList<String>> AATypeOptions;
    //COMETreeNode.assignments assigns each level an index in AATypeOptions.get(level), and thus an AA type
    //If -1, then no assignment yet

    int numMaxMut;//number of mutations allowed away from wtSeq (-1 means no cap)
    String wtSeq[];

    //information on states
    int numStates;//how many states there are
    //they have to have the same mutable residues & AA options,
    //though the residues involved may be otherwise different

    SearchProblemSuper[] mutableSearchProblems;//SearchProblems involved in COMETS search

    public SearchProblemSuper nonMutableSearchProblem;

    ArrayList<ArrayList<Integer>> mutable2StatePosNums;
    //mutable2StatePosNum.get(state) maps levels in this tree to flexible positions for state
    //(not necessarily an onto mapping)

    int stateNumPos[];

    int numSeqsReturned = 0;
    int stateGMECsForPruning = 0;//how many state GMECs have been calculated for nodes that are pruned    

    //HMN: To Help with Interaction E Heuristic
    int[] numMutPerStrand;
    //Maps the bound res num to the corresponding unbound emat
    HashMap<Integer, EnergyMatrix> boundResNumToUnboundEmat;
    //Maps the bound res num to the corresponding unbound res num
    HashMap<Integer, Integer> boundResNumToUnboundResNum;
    //Maps the bound res num to boolean that is true if res num is part of mutable
    //strand
    HashMap<Integer, Boolean> boundResNumToIsMutableStrand;
    //determines if two residues are on the same strand
    boolean[][] belongToSameStrand;

    public COMETSTreeSuper(int numTreeLevels, LME objFcn, LME[] constraints,
            ArrayList<ArrayList<String>> AATypeOptions, int numMaxMut, String[] wtSeq,
            int numStates, SearchProblemSuper[] stateSP, SearchProblemSuper aNonMutableSearchProblem,
            ArrayList<ArrayList<Integer>> mutable2StatePosNums) {

        this.numTreeLevels = numTreeLevels;
        this.objFcn = objFcn;
        this.constraints = constraints;
        this.AATypeOptions = AATypeOptions;
        this.numMaxMut = numMaxMut;
        this.wtSeq = wtSeq;
        this.numStates = numStates;
        this.mutableSearchProblems = stateSP;
        this.mutable2StatePosNums = mutable2StatePosNums;
        this.nonMutableSearchProblem = aNonMutableSearchProblem;
        mutableSearchProblems[0].emat.setConstTerm(0.0);
        mutableSearchProblems[1].emat.setConstTerm(0.0);
        nonMutableSearchProblem.emat.setConstTerm(0.0);
        
        stateNumPos = new int[numStates];
        for (int state = 0; state < numStates; state++) {
            stateNumPos[state] = stateSP[state].confSpaceSuper.numPos;
        }

        this.boundResNumToUnboundEmat = getBoundPosNumToUnboundEmat();
        this.boundResNumToUnboundResNum = getBoundPosNumToUnboundPosNum();
        this.boundResNumToIsMutableStrand = getBoundPosNumberToIsMutableStrand();
        this.belongToSameStrand = getSameStrandMatrix();
    }

    @Override
    public ArrayList<AStarNode> getChildren(AStarNode curNode) {
        COMETSNodeSuper seqNode = (COMETSNodeSuper) curNode;
        ArrayList<AStarNode> ans = new ArrayList<>();

        if (seqNode.isFullyDefined()) {
            seqNode.expandConfTree();
            ans.add(seqNode);
            return ans;
        } else {
            //expand next position...            
            int curAssignments[] = seqNode.getNodeAssignments();

            for (int splitPos = 0; splitPos < numTreeLevels; splitPos++) {
                if (curAssignments[splitPos] < 0) {//can split this level

                    for (int aa = 0; aa < AATypeOptions.get(splitPos).size(); aa++) {

                        int childAssignments[] = curAssignments.clone();
                        childAssignments[splitPos] = aa;

                        UpdatedPruningMatrixSuper childPruneMat[] = new UpdatedPruningMatrixSuper[numStates];
                        for (int state = 0; state < numStates; state++) {
                            childPruneMat[state] = doChildPruning(state, seqNode.pruneMat[state], splitPos, aa);
                        }

                        COMETSNodeSuper childNode = new COMETSNodeSuper(childAssignments, childPruneMat);

                        if (splitPos == numTreeLevels - 1) {//sequence now fully defined...make conf trees
                            makeSeqConfTrees(childNode);
                        }

                        childNode.setScore(boundLME(childNode, objFcn));

                        ans.add(childNode);
                    }

                    return ans;
                }
            }

            throw new RuntimeException("ERROR: Not splittable position found but sequence not fully defined...");
        }

    }

    private UpdatedPruningMatrixSuper doChildPruning(int state, PruningMatrix parentMat, int splitPos, int aa) {
        //Create an update to parentMat (without changing parentMat)
        //to reflect that splitPos has been assigned an amino-acid type

        String assignedAAType = AATypeOptions.get(splitPos).get(aa);

        UpdatedPruningMatrixSuper ans = new UpdatedPruningMatrixSuper(parentMat);
        int posAtState = mutable2StatePosNums.get(state).get(splitPos);

        //first, prune all other AA types at splitPos
        for (int rc : parentMat.unprunedRCsAtPos(posAtState)) {
            //HUNTER: TODO: AATYperPerRes should only be one residue for now
            //We should change this to allow for superRCs at the sequence search level
            String rcAAType = mutableSearchProblems[state].confSpaceSuper.posFlexSuper.get(posAtState).superRCs.get(rc).AATypePerRes.get(0);

            if (!rcAAType.equalsIgnoreCase(assignedAAType)) {
                ans.markAsPruned(new SuperRCTuple(posAtState, rc));
            }
        }

        //now do any consequent singles & pairs pruning
        int numUpdates = ans.countUpdates();
        int oldNumUpdates;

        PrunerSuper dee = new PrunerSuper(mutableSearchProblems[state], ans, true, Double.POSITIVE_INFINITY,
                0.0, false, mutableSearchProblems[state].useTupExpForSearch, false);
        //this is rigid, type-dependent pruning aiming for sequence GMECs
        //So Ew = Ival = 0

        do {//repeat as long as we're pruning things
            oldNumUpdates = numUpdates;
            dee.prune("GOLDSTEIN");
            dee.prune("GOLDSTEIN PAIRS FULL");
            numUpdates = ans.countUpdates();
        } while (numUpdates > oldNumUpdates);

        return ans;
    }

    private void makeSeqConfTrees(COMETSNodeSuper node) {
        //Given a node with a fully defined sequence, build its conformational search trees
        //for each state
        //If a state has no viable conformations, leave it null, with stateUB[state] = inf

        node.stateTrees = new ConfTreeSuper[numStates];
        node.stateUB = new double[numStates];

        for (int state = 0; state < numStates; state++) {

            //first make sure there are RCs available at each position
            boolean RCsAvailable = true;
            for (int pos = 0; pos < stateNumPos[state]; pos++) {
                if (node.pruneMat[state].unprunedRCsAtPos(pos).isEmpty()) {
                    RCsAvailable = false;
                    break;
                }
            }

            if (RCsAvailable) {
                node.stateTrees[state] = new ConfTreeSuper(mutableSearchProblems[state], node.pruneMat[state], false);

                AStarNode rootNode = node.stateTrees[state].rootNode();

                int blankConf[] = new int[stateNumPos[state]];//set up root node UB
                Arrays.fill(blankConf, -1);
                rootNode.UBConf = blankConf;
                updateUB(state, rootNode, node);
                node.stateUB[state] = rootNode.UB;

                node.stateTrees[state].initQueue(rootNode);//allocate queue and add root node
            } else {//no confs available for this state!
                node.stateTrees[state] = null;
                node.stateUB[state] = Double.POSITIVE_INFINITY;
            }
        }
    }

    @Override
    public AStarNode rootNode() {
        int[] conf = new int[numTreeLevels];
        Arrays.fill(conf, -1);//indicates sequence not assigned

        PruningMatrix[] pruneMats = new PruningMatrix[numStates];
        for (int state = 0; state < numStates; state++) {
            pruneMats[state] = mutableSearchProblems[state].pruneMat;
        }

        COMETSNodeSuper root = new COMETSNodeSuper(conf, pruneMats);
        root.setScore(boundLME(root, objFcn));
        return root;
    }

    @Override
    public boolean canPruneNode(AStarNode node) {
        //check if the node can be pruned based on the constraints
        //each constraint function must be <=0, so if any constraint function's lower bound
        //over sequences in seqNode is > 0,
        //we can prune seqNode
        COMETSNodeSuper seqNode = (COMETSNodeSuper) node;

        if (numMaxMut != -1) {
            //we have a cap on the number of mutations...prune if exceeded
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

            if (mutCount > numMaxMut)//No state GMECs calcd for this...pruning based only on sequence
            {
                return true;
            }
        }

        for (LME constr : constraints) {
            if (boundLME(seqNode, constr) > 0) {
                stateGMECsForPruning += countStateGMECs(seqNode);
                return true;
            }
        }

        return false;
    }

    @Override
    public boolean isFullyAssigned(AStarNode node) {
        //This checks if the node is returnable
        //So it must be fully processed (state GMECs found, not just fully defined sequence)

        if (!node.isFullyDefined())//sequence not fully defined
        {
            return false;
        }

        COMETSNodeSuper seqNode = (COMETSNodeSuper) node;

        for (int state = 0; state < numStates; state++) {
            if (seqNode.stateTrees[state] != null) {
                AStarNode bestNodeForState = seqNode.stateTrees[state].getQueue().peek();

                if (!bestNodeForState.isFullyDefined())//State GMEC calculation not done
                {
                    return false;
                }
            }
        }

        //if we get here, all state GMECs are calculated.  Node is fully processed
        return true;
    }

    @Override
    public int[] outputNode(AStarNode node) {
        //Let's print more info when outputting a node
        printBestSeqInfo((COMETSNodeSuper) node);
        return node.getNodeAssignments();
    }

    void printBestSeqInfo(COMETSNodeSuper seqNode) {
        //About to return the given fully assigned sequence from A*
        //provide information
        System.out.println("SeqTree: A* returning conformation; lower bound = " + seqNode.getScore() + " nodes expanded: " + numExpanded);

        numSeqsReturned++;

        System.out.print("Sequence: ");

        for (int level = 0; level < numTreeLevels; level++) {
            System.out.print(AATypeOptions.get(level).get(seqNode.getNodeAssignments()[level]) + " ");
        }
        System.out.println();

        //provide state GMECs, specified as rotamers (AA types all the same of course)
        for (int state = 0; state < numStates; state++) {
            System.out.print("State " + state);

            if (seqNode.stateTrees[state] == null) {
                System.out.println(" has an unavoidable clash.");
            } else {
                System.out.print(" superRCs: ");
                int conf[] = seqNode.stateTrees[state].getQueue().peek().getNodeAssignments();
                for (int pos = 0; pos < stateNumPos[state]; pos++) {
                    System.out.print(conf[pos] + " ");
                }

                System.out.println("Energy: "
                        + getEnergyMatrix(state).confE(conf));
            }
        }
        int[] confBound = seqNode.stateTrees[0].getQueue().peek().getNodeAssignments();
        int[] confUnbound = seqNode.stateTrees[1].getQueue().peek().getNodeAssignments();
        double LMEscore = getEnergyMatrix(0).confE(confBound) - getEnergyMatrix(1).confE(confUnbound);
        System.out.println("LME: " + LMEscore);

        int numSeqDefNodes = 0;
        for (AStarNode node : getQueue()) {
            if (node.isFullyDefined()) {
                numSeqDefNodes++;
            }
        }

        System.out.println();
        System.out.println(numExpanded + " expanded; " + getQueue().size() + " nodes in tree, of which "
                + numSeqDefNodes + " are fully defined; "
                + numPruned + " pruned.");

        int stateGMECsRet = numSeqsReturned * numStates;
        int stateGMECsInTree = countGMECsInTree();
        int totGMECsCalcd = stateGMECsRet + stateGMECsInTree + stateGMECsForPruning;
        System.out.println(totGMECsCalcd + " state GMECs calculated: " + stateGMECsRet + " returned, " + stateGMECsInTree
                + " in tree, " + stateGMECsForPruning + " for pruned sequences.");

    }

    int countGMECsInTree() {
        //count how many state GMECs have been calculated and are in nodes in the tree--
        //that is, how many ConfTrees at sequence nodes have been expanded to the points
        //that their lowest-bound node is fully expanded (and thus is their GMEC)
        //for comparison, regular A* needs to calculate the GMEC for every (non-clashing)
        //(sequence,state) pair
        //Note this only counts state GMECs currently in the tree
        int count = 0;

        for (AStarNode node : getQueue()) {
            count += countStateGMECs((COMETSNodeSuper) node);
        }

        return count;
    }

    int countStateGMECs(COMETSNodeSuper seqNode) {
        //how many states for this node have GMECs calculated?
        int count = 0;

        if (seqNode.stateTrees != null) {//fully defined sequence, so there are state trees
            for (ConfTreeSuper ct : seqNode.stateTrees) {
                if (ct != null) {
                    if (ct.getQueue().peek() != null) {
                        AStarNode bestNode = ct.getQueue().peek();
                        if (bestNode.isFullyDefined()) {
                            count++;
                        }
                    }
                }
            }
        }

        return count;
    }

    /**
     *
     * @param seqNode the sequence node in the A* tree
     * @param func the LME describing the objective function
     * @return lower bound on the objective function
     */
    private double boundLME(COMETSNodeSuper seqNode, LME func) {
        //Lower-bound func over the sequence space defined by this node
        if (seqNode.isFullyDefined())//fully-defined sequence
        {
            double originalBound = calcLBConfTrees(seqNode, func);
            //double exactBound = calcLBConfTreeMPLP(seqNode);
            return originalBound;
        } else {

            double newBound = calcLBPartialSeqImproved(seqNode);
            //double maxInterfaceBound = maxInterfaceBound(seqNode);
            //double newBoundV2 = calcLBPartialSeqImprovedVersion2(seqNode);
            double maxInterfaceBoundImproved = maxInterfaceBoundWithProtein(seqNode);
            double originalBound = calcLBPartialSeq(seqNode, func);
            printSequence(getSequence(seqNode));
            System.out.println("New bound: " + newBound);
//            System.out.println("Max Interface Improved: " + maxInterfaceBoundImproved);
            System.out.println("Original bound: " + originalBound);
            System.out.println();
            //return Math.max(Math.max(maxInterfaceBoundImproved, newBound), originalBound);
            return newBound;
        }
    }

    private int[] AAOptions(int pos, int assignment) {
        //which of the AA option indices for this mutable position are allowed 
        //given this assignment in nodeAssignments?

        if (assignment == -1) {//all options allowed
            int numOptions = AATypeOptions.get(pos).size();
            int[] ans = new int[numOptions];
            for (int option = 0; option < numOptions; option++) {
                ans[option] = option;
            }

            return ans;
        } else//just the assigned option
        {
            return new int[]{assignment};
        }
    }

    private double calcLBPartialSeq(COMETSNodeSuper seqNode, LME func) {

        int partialSeq[] = seqNode.getNodeAssignments();

        double ans = func.constTerm;
        //first terms for mutable residues
        for (int i = 0; i < numTreeLevels; i++) {
            double resE = Double.POSITIVE_INFINITY;

            for (int curAA : AAOptions(i, partialSeq[i])) {
                double AAE = 0;

                //get contributions to residue energy for this AA type from all states
                for (int state = 0; state < numStates; state++) {
                    if (func.coeffs[state] != 0) {

                        int statePosNum = mutable2StatePosNums.get(state).get(i);//residue number i converted to this state's flexible position numbering
                        boolean minForState = (func.coeffs[state] > 0);//minimizing (instead of maximizing) energy for this state

                        double stateAAE = Double.POSITIVE_INFINITY;

                        PruningMatrix pruneMat = seqNode.pruneMat[state];
                        EnergyMatrix eMatrix = getEnergyMatrix(state);

                        ArrayList<Integer> rotList = unprunedRCsAtAA(state, pruneMat,
                                i, statePosNum, curAA);

                        if (!minForState) {
                            stateAAE = Double.NEGATIVE_INFINITY;
                        }

                        for (int rot : rotList) {

                            double rotE = eMatrix.getOneBody(statePosNum, rot);

                            for (int pos2 = 0; pos2 < stateNumPos[state]; pos2++) {//all non-mut; seq only if < this one
                                if ((!mutable2StatePosNums.get(state).contains(pos2)) || pos2 < statePosNum) {

                                    double bestInteraction = Double.POSITIVE_INFINITY;
                                    ArrayList<Integer> rotList2 = pruneMat.unprunedRCsAtPos(pos2);

                                    if (!minForState) {
                                        bestInteraction = Double.NEGATIVE_INFINITY;
                                    }

                                    for (int rot2 : rotList2) {
                                        //rot2 known to be unpruned
                                        if (!pruneMat.getPairwise(statePosNum, rot, pos2, rot2)) {

                                            double pairwiseE = eMatrix.getPairwise(statePosNum, rot, pos2, rot2);
                                            pairwiseE += higherOrderContrib(state, pruneMat, statePosNum, rot, pos2, rot2, minForState);

                                            if (minForState) {
                                                bestInteraction = Math.min(bestInteraction, pairwiseE);
                                            } else {
                                                bestInteraction = Math.max(bestInteraction, pairwiseE);
                                            }
                                        }
                                    }

                                    rotE += bestInteraction;
                                }
                            }

                            if (minForState) {
                                stateAAE = Math.min(stateAAE, rotE);
                            } else {
                                stateAAE = Math.max(stateAAE, rotE);
                            }
                        }

                        if (Double.isInfinite(stateAAE)) {
                            //this will occur if the state is impossible (all confs pruned)
                            if (func.coeffs[state] > 0) {
                                AAE = Double.POSITIVE_INFINITY;
                            } else if (func.coeffs[state] < 0 && AAE != Double.POSITIVE_INFINITY) //if a "positive-design" (coeff>0) state is impossible, return +inf overall
                            {
                                AAE = Double.NEGATIVE_INFINITY;
                            }
                            //else AAE unchanged, since func doesn't involve this state
                        } else {
                            AAE += func.coeffs[state] * stateAAE;
                        }
                    }
                }

                resE = Math.min(resE, AAE);
            }

            ans += resE;
        }

        //now we bound the energy for the other residues for each of the states
        //(internal energy for that set of residues)
        for (int state = 0; state < numStates; state++) {

            if (func.coeffs[state] != 0) {
                ans += func.coeffs[state] * getEnergyMatrix(state).getConstTerm();

                double nonMutBound = boundStateNonMutE(state, seqNode, func.coeffs[state] > 0);
                //handle this like AAE above
                if (Double.isInfinite(nonMutBound)) {
                    //this will occur if the state is impossible (all confs pruned)
                    if (func.coeffs[state] > 0) {
                        return Double.POSITIVE_INFINITY;
                    } else if (func.coeffs[state] < 0) {
                        ans = Double.NEGATIVE_INFINITY;
                    }
                    //if a "positive-design" (coeff>0) state is impossible, return +inf overall still
                    //else ans unchanged, since func doesn't involve this state
                } else {
                    ans += func.coeffs[state] * nonMutBound;
                }
            }
        }

        return ans;
    }

    double higherOrderContrib(int state, PruningMatrix pruneMat,
            int pos1, int rc1, int pos2, int rc2, boolean minForState) {
        //higher-order contribution for a given RC pair in a given state, 
        //when scoring a partial conf

        EnergyMatrix emat = getEnergyMatrix(state);
        HigherTupleFinder<Double> htf = emat.getHigherOrderTerms(pos1, rc1, pos2, rc2);

        if (htf == null) {
            return 0;//no higher-order interactions
        } else {
            SuperRCTuple curPair = new SuperRCTuple(pos1, rc1, pos2, rc2);
            return higherOrderContrib(state, pruneMat, htf, curPair, minForState);
        }
    }

    private boolean posComesBefore(int pos1, int pos2, int state) {
        //Does pos1 come "before" pos2?  (Ordering: non-mut pos, then mut pos, in ascending order)
        //These are flexible positions for the given state
        boolean isMut1 = mutable2StatePosNums.get(state).contains(pos1);
        boolean isMut2 = mutable2StatePosNums.get(state).contains(pos2);

        if (isMut1 && !isMut2) {
            return false;
        }
        if (isMut2 && !isMut1) {
            return true;
        }

        //ok they are in the same category (mut or not)
        return pos1 < pos2;
    }

    double higherOrderContrib(int state, PruningMatrix pruneMat, HigherTupleFinder<Double> htf,
            SuperRCTuple startingTuple, boolean minForState) {
        //recursive function to get bound on higher-than-pairwise terms
        //this is the contribution to the bound due to higher-order interactions
        //of the RC tuple startingTuple (corresponding to htf)

        double contrib = 0;

        //to avoid double-counting, we are just counting interactions of starting tuple
        //with residues before the "earliest" one (startingLevel) in startingTuple
        //"earliest" means lowest-numbered, except non-mutating res come before mutating
        int startingLevel = startingTuple.pos.get(startingTuple.pos.size() - 1);

        for (int iPos : htf.getInteractingPos()) {//position has higher-order interaction with tup
            if (posComesBefore(iPos, startingLevel, state)) {//interaction in right order
                //(want to avoid double-counting)

                double levelBestE = Double.POSITIVE_INFINITY;//best value of contribution
                //from tup-iPos interaction

                if (!minForState) {
                    levelBestE = Double.NEGATIVE_INFINITY;
                }

                ArrayList<Integer> allowedRCs = pruneMat.unprunedRCsAtPos(iPos);

                for (int rc : allowedRCs) {

                    SuperRCTuple augTuple = startingTuple.addRC(iPos, rc);

                    if (!pruneMat.isPruned(augTuple)) {

                        double interactionE = htf.getInteraction(iPos, rc);

                        //see if need to go up to highers order again...
                        HigherTupleFinder htf2 = htf.getHigherInteractions(iPos, rc);
                        if (htf2 != null) {
                            interactionE += higherOrderContrib(state, pruneMat, htf2, augTuple, minForState);
                        }

                        //besides that only residues in definedTuple or levels below level2
                        if (minForState) {
                            levelBestE = Math.min(levelBestE, interactionE);
                        } else {
                            levelBestE = Math.max(levelBestE, interactionE);
                        }
                    }
                }

                contrib += levelBestE;//add up contributions from different interacting positions iPos
            }
        }

        return contrib;
    }

    private ArrayList<Integer> unprunedRCsAtAA(int state, PruningMatrix pruneMat,
            int mutablePos, int statePosNum, int curAA) {
        //List the RCs of the given position in the given state
        //that come with the indicated AA type (indexed in AATypeOptions)

        ArrayList<Integer> unprunedRCs = pruneMat.unprunedRCsAtPos(statePosNum);
        ArrayList<Integer> ans = new ArrayList<>();

        String curAAType = AATypeOptions.get(mutablePos).get(curAA);

        for (int rc : unprunedRCs) {
            //HMN: TODO Change to allow for more than one AA type per superRCs.
            //Currently this does not support merging
            String rcAAType = mutableSearchProblems[state].confSpaceSuper.posFlexSuper.get(statePosNum).superRCs.get(rc).AATypePerRes.get(0);

            if (rcAAType.equalsIgnoreCase(curAAType))//right type
            {
                ans.add(rc);
            }
        }

        return ans;
    }


    private double calcLBConfTrees(COMETSNodeSuper seqNode, LME func) {
        //here the sequence is fully defined
        //so we can bound func solely based on lower and upper bounds (depending on func coefficient sign)
        //of the GMECs for each state, which can be derived from the front node of each state's ConfTree
        double ans = func.constTerm;
        for (int state = 0; state < numStates; state++) {
            if (func.coeffs[state] > 0) {//need lower bound

                if (seqNode.stateTrees[state] == null)//state and sequence impossible
                {
                    return Double.POSITIVE_INFINITY;
                }

                ans += func.coeffs[state] * seqNode.stateTrees[state].getQueue().peek().getScore();
            } else if (func.coeffs[state] < 0) {//need upper bound

                ConfTreeSuper curTree = seqNode.stateTrees[state];

                if (curTree == null) {//state and sequence impossible
                    //bound would be +infinity
                    ans = Double.NEGATIVE_INFINITY;
                    continue;
                }

                //make sure stateUB is updated, at least based on the current best node in this state's tree
                PriorityQueue<AStarNode> curExpansion = curTree.getQueue();

                AStarNode curNode = curExpansion.peek();

                if (Double.isNaN(curNode.UB)) {//haven't calculated UB, so calculate and update stateUB
                    while (!updateUB(state, curNode, seqNode)) {
                        //if UB not calc'd yet, curNode.UBConf is from curNode's parent
                        //if updateUB fails then the node has an inevitable clash...remove it from the tree
                        curExpansion.poll();
                        if (curExpansion.isEmpty()) {//no more nodes, so no available states!
                            //bound would be +infinity
                            seqNode.stateUB[state] = Double.POSITIVE_INFINITY;
                            ans = Double.NEGATIVE_INFINITY;
                            curNode = null;
                            break;
                        }

                        curNode = curExpansion.peek();
                        if (!Double.isNaN(curNode.UB)) {
                            break;
                        }
                    }
                }

                if (curNode == null) {//no nodes left for this state
                    ans = Double.NEGATIVE_INFINITY;
                    continue;
                }

                seqNode.stateUB[state] = Math.min(seqNode.stateUB[state], curNode.UB);

                ans += func.coeffs[state] * seqNode.stateUB[state];
            }

            //shell-shell energies may differ between states!
            //THIS DOUBLE-COUNTS BC SCORES ALREADY INCLUDE CONST TERM
            //ans += func.coeffs[state]*getEnergyMatrix(state).getConstTerm();
        }

        return ans;
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

    private double boundStateNonMutE(int state, COMETSNodeSuper seqNode, boolean minForState) {
        //get a quick lower or upper bound (as indicated) for the energy of the given state's
        //non-mutable residues (their intra+shell energies + pairwise between them)
        //use pruning information from seqNode
        double ans = 0;

        PruningMatrix pruneMat = seqNode.pruneMat[state];
        EnergyMatrix eMatrix = getEnergyMatrix(state);

        for (int pos = 0; pos < stateNumPos[state]; pos++) {
            if ((!mutable2StatePosNums.get(state).contains(pos))) {

                double resE = Double.POSITIVE_INFINITY;

                ArrayList<Integer> rotList = pruneMat.unprunedRCsAtPos(pos);

                if (!minForState) {
                    resE = Double.NEGATIVE_INFINITY;
                }

                for (int rot : rotList) {
                    //make sure rot isn't pruned
                    if (!pruneMat.getOneBody(pos, rot)) {

                        double rotE = eMatrix.getOneBody(pos, rot);

                        for (int pos2 = 0; pos2 < stateNumPos[state]; pos2++) {//all non-mut; seq only if < this one
                            if ((!mutable2StatePosNums.get(state).contains(pos2)) && pos2 < pos) {

                                double bestInteraction = Double.POSITIVE_INFINITY;
                                ArrayList<Integer> rotList2 = pruneMat.unprunedRCsAtPos(pos2);
                                if (!minForState) {
                                    bestInteraction = Double.NEGATIVE_INFINITY;
                                }

                                for (int rot2 : rotList2) {
                                    if (!pruneMat.getPairwise(pos, rot, pos2, rot2)) {
                                        double pairwiseE = eMatrix.getPairwise(pos, rot, pos2, rot2);

                                        pairwiseE += higherOrderContrib(state, pruneMat, pos, rot, pos2, rot2, minForState);

                                        if (minForState) {
                                            bestInteraction = Math.min(bestInteraction, pairwiseE);
                                        } else {
                                            bestInteraction = Math.max(bestInteraction, pairwiseE);
                                        }
                                    }
                                }

                                rotE += bestInteraction;
                            }
                        }

                        if (minForState) {
                            resE = Math.min(resE, rotE);
                        } else {
                            resE = Math.max(resE, rotE);
                        }
                    }
                }

                ans += resE;
            }
        }

        return ans;
    }

    boolean updateUB(int state, AStarNode expNode, COMETSNodeSuper seqNode) {
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
                posOptions = seqNode.pruneMat[state].unprunedRCsAtPos(pos);
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

    //For flexible position in bound matrix (0,1,2,...,numRes-1) we map to the corresponding
    //unbound energy matrix
    public HashMap<Integer, EnergyMatrix> getBoundPosNumToUnboundEmat() {
        SearchProblemSuper boundState = this.mutableSearchProblems[0];
        //Get res number from each flexible position in the bound state
        List<String> resNumsBound = boundState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound mutable state
        SearchProblemSuper unBoundMutableState = this.mutableSearchProblems[1];
        List<String> resNumsUnboundMutable = unBoundMutableState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound non-mutable state
        SearchProblemSuper unBoundNonMutableState = this.nonMutableSearchProblem;
        List<String> resNumsUnboundNonMutable = unBoundNonMutableState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
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

    //For flexible position in bound matrix (0,1,2,...,numRes-1) we map to the corresponding
    //position number in the unbound matrix
    public HashMap<Integer, Integer> getBoundPosNumToUnboundPosNum() {
        SearchProblemSuper boundState = this.mutableSearchProblems[0];
        //Get res number from each flexible position in the bound state
        List<String> resNumsBound = boundState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound mutable state
        SearchProblemSuper unBoundMutableState = this.mutableSearchProblems[1];
        List<String> resNumsUnboundMutable = unBoundMutableState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound non-mutable state
        SearchProblemSuper unBoundNonMutableState = this.nonMutableSearchProblem;
        List<String> resNumsUnboundNonMutable = unBoundNonMutableState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
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

    //For flexible position in bound matrix (0,1,2,...,numRes-1) we map to the corresponding
    //strand number
    public boolean[][] getSameStrandMatrix() {
        SearchProblemSuper boundState = this.mutableSearchProblems[0];
        //Get res number from each flexible position in the bound state
        List<String> resNumsBound = boundState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toList());

        //Get res number for each flexible position in the unbound mutable state
        SearchProblemSuper unBoundMutableState = this.mutableSearchProblems[1];
        List<String> resNumsUnboundMutable = unBoundMutableState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toList());

        //Get res number for each flexible position in the unbound non-mutable state
        SearchProblemSuper unBoundNonMutableState = this.nonMutableSearchProblem;
        List<String> resNumsUnboundNonMutable = unBoundNonMutableState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
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

    //For flexible position in bound matrix (0,1,2,...,numRes-1) we map to the boolean
    //that is true if the res is part of boolean strand and false otherwise
    public HashMap<Integer, Boolean> getBoundPosNumberToIsMutableStrand() {
        SearchProblemSuper boundState = this.mutableSearchProblems[0];
        //Get res number from each flexible position in the bound state
        List<String> resNumsBound = boundState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound mutable state
        SearchProblemSuper unBoundMutableState = this.mutableSearchProblems[1];
        List<String> resNumsUnboundMutable = unBoundMutableState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound non-mutable state
        SearchProblemSuper unBoundNonMutableState = this.nonMutableSearchProblem;
        List<String> resNumsUnboundNonMutable = unBoundNonMutableState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
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
     * calcLBPartialSeqImproved: Computes a lower bound on a multi state energy
     * of a partial sequence assignments for a PROTEIN:LIGAND interaction. Our
     * bound consists of (in BOLTZMANN weighted terms):
     * MAX(P,LA,P:LA)/(MAX(P)*MAX(LA)) *
     * MAX_S((MAX(P:LU_s)*MAX(LA:LU_s))/(MIN(LA:LU_s))) where P is the target
     * protein whose sequence is known, LA are the ligand assigned residues
     * whose sequence has been defined in seqNode, and LU_s are the ligand
     * unassigned residues In ENERGIES, our bound is: GMinEC(P,LA,P:LA) -
     * GMinEC(P) - GMinEC(LA) + MIN_S(GMinEC(P:LU_s) + GMinEC(LA:LU_s) -
     * GMaxEC(LA:LU_s)) where GMinEC is the energy of the global minimum energy
     * conformation and GMaxEC is the energy of the global maximum energy
     * conformation
     *
     * @param seqNode
     * @param boundResNumToUnboundEmat
     * @param boundResNumToUnboundResNum
     * @param boundResNumToIsMutableStrand
     * @param belongToSameStrand
     * @return
     */
    private double calcLBPartialSeqImproved(COMETSNodeSuper seqNode) {
        SearchProblemSuper boundSP = mutableSearchProblems[0];
        SearchProblemSuper ligandSP = mutableSearchProblems[1];
        SearchProblemSuper proteinSP = nonMutableSearchProblem;

        // First compute GMinEC(P,LA,P:LA). Here an lower bound can be used, but ideally it should be computed exactly
        // This term involves the bound state
        double gminec_p_la_pla = 0;
        ArrayList<Integer> subsetOfPositions_p_la_pla = getLigandAssignedPosNums(seqNode, true);
        subsetOfPositions_p_la_pla.addAll(getProteinPosNums(true));
        Collections.sort(subsetOfPositions_p_la_pla);

        SearchProblemSuper searchProblem_p_la_pla = boundSP.getSubsetSearchProblem(subsetOfPositions_p_la_pla);
        updateSubsetPruneMat(searchProblem_p_la_pla.pruneMat, seqNode.pruneMat[0], subsetOfPositions_p_la_pla);
        ConfTreeSuper confTree_p_la_pla = new ConfTreeSuper(searchProblem_p_la_pla);
        //ConfTreeSuper confTree_p_la_pla = getPartialConfTree_InternalTerm(boundSP, subsetOfPositions_p_la_pla, seqNode.pruneMat[0]);
        gminec_p_la_pla = confTree_p_la_pla.energyNextConf();

        // GMinEC(P) can be precomputed because it is a constant for the system or computed here. 
        // This term involves the unbound protein
        double gminec_p = 0;
        ArrayList<Integer> subsetOfPositions_p = getProteinPosNums(false);
        // Technically, we could just create a new ConfTreeSuper() directly, since we are using the entire protein
        // However, we also want to make sure we don't have any shell-shell energies in this calculation, so 
        //  I will just use the getPartialConfTree() method
        SearchProblemSuper searchProblem_p = proteinSP.getSubsetSearchProblem(subsetOfPositions_p);
        ConfTreeSuper confTree_p = new ConfTreeSuper(searchProblem_p);
        gminec_p = confTree_p.energyNextConf();

        // Now compute GMinEC(LA). This has to be computed exactly (or through a lower bound)
        // This term involves the unbound ligand 
        double gminec_la = 0;
        ArrayList<Integer> subsetofPositions_la = getLigandAssignedPosNums(seqNode, false);

        SearchProblemSuper searchProblem_la = ligandSP.getSubsetSearchProblem(subsetofPositions_la);
        updateSubsetPruneMat(searchProblem_la.pruneMat, seqNode.pruneMat[1], subsetofPositions_la);
        ConfTreeSuper confTree_la = new ConfTreeSuper(searchProblem_la);
        gminec_la = confTree_la.energyNextConf();

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Finally, compute MIN_S(GMinEC(P:LU_s) + GMinEC(LA:LU_s) - GMaxEC(LA:LU_s)).        
        // The following section should be "modular" because there are two ways to do this. One way is using a greedy algorithm to compute it
        // exactly. We have not developed it yet. For now let's compute it t as follows :
        //   MIN_S(GMinEC(P:LU_s) + GMinEC(LA:LU_s) - GMaxEC(LA:LU_s)) <= MIN_S(GMinEC(P:LU_s) + GMinEC(LA:LU_s)) - MAX_S(GMaxEC(LA:LU_s)
        // Thus, first compute: MIN_S(GMinEC(P:LU_s) + GMinEC(LA:LU_s)), which can be easily computed by computing a min gmec over: 
        //		GMinEC(LA:LU_s, P:LU_s), including all rotamers for all amino acids defined in s.
        //This term involves the bound state
        double gminec_lalus_plus = 0;
        ArrayList<Integer> ligandAssignedBound = getLigandAssignedPosNums(seqNode, true);
        ArrayList<Integer> ligandUnassignedBound = getLigandUnassignedPosNums(seqNode, true);
        ArrayList<Integer> proteinBound = getProteinPosNums(true);
        ArrayList<Integer> subsetOfPositions_lalus_plus = new ArrayList<>();
        subsetOfPositions_lalus_plus.addAll(ligandAssignedBound);
        subsetOfPositions_lalus_plus.addAll(ligandUnassignedBound);
        subsetOfPositions_lalus_plus.addAll(proteinBound);
        Collections.sort(subsetOfPositions_lalus_plus);
        // We create an interaction graph for LA:LU_s and for P:LU_s and then add them into one complete interaction graph
        boolean[][] interactionGraph__lalus_plus = addInteractionGraphs(createInteractionGraph(subsetOfPositions_lalus_plus, ligandAssignedBound, ligandUnassignedBound), createInteractionGraph(subsetOfPositions_lalus_plus, proteinBound, ligandUnassignedBound));

        SearchProblemSuper searchProblem_lalus_plus = boundSP.getSubsetSearchProblem(subsetOfPositions_lalus_plus);
        updateSubsetPruneMat(searchProblem_lalus_plus.pruneMat, seqNode.pruneMat[0], subsetOfPositions_lalus_plus);
        updateEnergyMatrix_CrossTerms(searchProblem_lalus_plus, subsetOfPositions_lalus_plus, interactionGraph__lalus_plus, false);
        addCrossTermInternalEnergies(searchProblem_lalus_plus.emat, boundSP.emat, ligandSP.emat, getBoundPosNumToUnboundPosNum(), subsetOfPositions_lalus_plus, ligandUnassignedBound);
        ConfTreeSuper confTree_lalus_plus_2 = new ConfTreeSuper(searchProblem_lalus_plus);
        gminec_lalus_plus = confTree_lalus_plus_2.energyNextConf();

        // Then compute the maximum MAX_S(GMaxEC(LA:LU_s), which can be computed by either negating all the energies in the matrix or something similar.
        // This involves an unbound state
        double gmaxec_lalus = 0;
        ArrayList<Integer> ligandAssignedUnbound = getLigandAssignedPosNums(seqNode, false);
        ArrayList<Integer> ligandUnassignedUnbound = getLigandUnassignedPosNums(seqNode, false);
        ArrayList<Integer> subsetOfPositions_lalus = new ArrayList<>();
        subsetOfPositions_lalus.addAll(ligandAssignedUnbound);
        subsetOfPositions_lalus.addAll(ligandUnassignedUnbound);
        Collections.sort(subsetOfPositions_lalus);

        boolean[][] interactionGraph_lalus = createInteractionGraph(subsetOfPositions_lalus, ligandAssignedUnbound, ligandUnassignedUnbound);
        //create confTree and negate energies

        SearchProblemSuper searchProblem = ligandSP.getSubsetSearchProblem(subsetOfPositions_lalus);
        updateSubsetPruneMat(searchProblem.pruneMat, seqNode.pruneMat[1], subsetOfPositions_lalus);
        updateEnergyMatrix_CrossTerms(searchProblem, subsetOfPositions_lalus, interactionGraph_lalus, true);
        ConfTreeSuper confTree_lalus_2 = new ConfTreeSuper(searchProblem);
        gmaxec_lalus = -confTree_lalus_2.energyNextConf();

        double lowerBound = gminec_p_la_pla - gminec_p - gminec_la + gminec_lalus_plus - gmaxec_lalus;
        lowerBound += this.mutableSearchProblems[0].emat.getConstTerm() - this.mutableSearchProblems[1].emat.getConstTerm() - this.nonMutableSearchProblem.emat.getConstTerm();

        return lowerBound;

    }

    private double maxInterfaceBoundWithProtein(COMETSNodeSuper seqNode) {
        SearchProblemSuper boundSP = mutableSearchProblems[0];
        SearchProblemSuper ligandSP = mutableSearchProblems[1];
        SearchProblemSuper proteinSP = nonMutableSearchProblem;

        ArrayList<Integer> proteinPosNums = getProteinPosNums(true);
        ArrayList<Integer> ligandPosNums = getLigandAssignedPosNums(seqNode, true);
        ligandPosNums.addAll(getLigandUnassignedPosNums(seqNode, true));
        Collections.sort(ligandPosNums);

        ArrayList<Integer> allPos = new ArrayList<>();
        allPos.addAll(ligandPosNums);
        allPos.addAll(proteinPosNums);
        Collections.sort(allPos);

        boolean[][] interactionGraph = createInteractionGraph(allPos, ligandPosNums, proteinPosNums);

        SearchProblemSuper searchProblem_interface = boundSP.getSubsetSearchProblem(allPos);
        updateSubsetPruneMat(searchProblem_interface.pruneMat, seqNode.pruneMat[0], allPos);
        subtractInternalEnergies(searchProblem_interface.emat, boundSP.emat, allPos, ligandPosNums);
        subtractPairwiseEnergies(searchProblem_interface.emat, boundSP.emat, allPos, ligandPosNums);
        addCrossTermInternalEnergies(searchProblem_interface.emat, boundSP.emat, ligandSP.emat, boundResNumToUnboundResNum, allPos, ligandPosNums);
        ConfTreeSuper confTree_interface = new ConfTreeSuper(searchProblem_interface);
        double gmecInterface = confTree_interface.energyNextConf();

        double constTerms = boundSP.emat.getConstTerm() - ligandSP.emat.getConstTerm() + objFcn.constTerm;
        return gmecInterface + constTerms;
    }

    /**
     * We call cross-term internal energy (internalEnergyBound -
     * internalEnergyUnbound) This method add cross-term internal energies to
     * the confTree energy matrix
     *
     * @param confTree the confTree whose energy matrix we are going to change
     * @param origBoundEmat the original bound energy matrix (from search
     * problem)
     * @param origUnboundEmat the original unbound energy matrix (from search
     * problem)
     * @param boundToUnboundPosNum mapping from bound pos number to unbound pos
     * number to index into unbound energy matrix
     * @param allPositions all the positions in this confTree
     * @param posToAddInternalE the positions in this confTree for which we are
     * adding the cross internal energy
     */
    private void addCrossTermInternalEnergies(EnergyMatrix currentEmat, EnergyMatrix origBoundEmat, EnergyMatrix origUnboundEmat, HashMap<Integer, Integer> boundToUnboundPosNum, ArrayList<Integer> allPositions, ArrayList<Integer> posToAddInternalE) {
        for (int pos : posToAddInternalE) {
            //pos is the index in the old emat
            int unboundPos = boundToUnboundPosNum.get(pos);
            //get the index of this position in the confTree emat
            int indexCurrentEmat = allPositions.indexOf(pos);

            for (int rot = 0; rot < currentEmat.oneBody.get(indexCurrentEmat).size(); rot++) {
                double internalEBound = origBoundEmat.getOneBody(pos, rot);
                double internalEUnbound = origUnboundEmat.getOneBody(unboundPos, rot);
                double newInternalE = internalEBound - internalEUnbound;
                double currentIntraE = currentEmat.getOneBody(indexCurrentEmat, rot);
                currentEmat.setOneBody(pos, rot, newInternalE + currentIntraE);
            }
        }
    }

    private void addInternalEnergies(EnergyMatrix currentEmat, EnergyMatrix origEmat, ArrayList<Integer> allPositions, ArrayList<Integer> posToAddInternal) {
        addSubtractInternalEnergies(currentEmat, origEmat, allPositions, posToAddInternal, true);
    }

    private void subtractInternalEnergies(EnergyMatrix currentEmat, EnergyMatrix origEmat, ArrayList<Integer> allPositions, ArrayList<Integer> posToAddInternal) {
        addSubtractInternalEnergies(currentEmat, origEmat, allPositions, posToAddInternal, false);
    }

    /**
     * This method will either add or subtract internal energies in a confTree
     * energy matrix
     *
     * @param currentEmat the energy matrix we are updating
     * @param origEmat the energy matrix from which we are getting energies to
     * add/subtract
     * @param allPositions the positions in the confTree emat (each int is the
     * original posNum from the original search problem)
     * @param posToAddInternal the positions that we want to add internal
     * energies for
     * @param add should we add or subtract energies
     */
    private void addSubtractInternalEnergies(EnergyMatrix currentEmat, EnergyMatrix origEmat, ArrayList<Integer> allPositions, ArrayList<Integer> posToAddInternal, boolean add) {
        for (int pos : posToAddInternal) {
            int indexConfTreeEmat = allPositions.indexOf(pos);
            for (int rot = 0; rot < currentEmat.oneBody.get(indexConfTreeEmat).size(); rot++) {
                double currentE = currentEmat.getOneBody(indexConfTreeEmat, rot);
                double toAddOrSubtract = origEmat.getOneBody(pos, rot);
                if (add) {
                    currentEmat.setOneBody(indexConfTreeEmat, rot, currentE + toAddOrSubtract);
                } //else we subtract
                else {
                    currentEmat.setOneBody(indexConfTreeEmat, rot, currentE - toAddOrSubtract);
                }
            }
        }
    }

    private void addPairwiseEnergies(EnergyMatrix currentEmat, EnergyMatrix origEmat, ArrayList<Integer> allPositions, ArrayList<Integer> posToAddPairwise) {
        addSubtractPairwiseEnergies(currentEmat, origEmat, allPositions, posToAddPairwise, true);
    }

    private void subtractPairwiseEnergies(EnergyMatrix currentEmat, EnergyMatrix origEmat, ArrayList<Integer> allPositions, ArrayList<Integer> posToAddPairwise) {
        addSubtractPairwiseEnergies(currentEmat, origEmat, allPositions, posToAddPairwise, false);
    }

    /**
     * This method will either add or subtract pairwise energies in a confTree
     * energy matrix
     *
     * @param confTree the confTree whose energy matrix we are adapting
     * @param origEmat the energy matrix from which we are getting energies to
     * add/subtract
     * @param allPositions the positions in the confTree emat (each int is the
     * original posNum from the original search problem)
     * @param posToAddInternal the positions that we want to add pairwise
     * energies for
     * @param add should we add or subtract energies
     */
    private void addSubtractPairwiseEnergies(EnergyMatrix currentEmat, EnergyMatrix origEmat, ArrayList<Integer> allPositions, ArrayList<Integer> posToAddPairwise, boolean add) {
        for (int i = 0; i < posToAddPairwise.size(); i++) {
            for (int j = i + 1; j < posToAddPairwise.size(); j++) {
                int indexConfTreeEmat_I = allPositions.indexOf(posToAddPairwise.get(i));
                int indexConfTreeEmat_J = allPositions.indexOf(posToAddPairwise.get(j));
                for (int rotI = 0; rotI < currentEmat.oneBody.get(indexConfTreeEmat_I).size(); rotI++) {
                    for (int rotJ = 0; rotJ < currentEmat.oneBody.get(indexConfTreeEmat_J).size(); rotJ++) {
                        double currentE = currentEmat.getPairwise(indexConfTreeEmat_I, rotI, indexConfTreeEmat_J, rotJ);
                        double toAddSubtract = origEmat.getPairwise(posToAddPairwise.get(i), rotI, posToAddPairwise.get(j), rotJ);
                        if (add) {
                            currentEmat.setPairwise(indexConfTreeEmat_I, rotI, indexConfTreeEmat_J, rotJ, currentE + toAddSubtract);
                        } else {
                            currentEmat.setPairwise(indexConfTreeEmat_I, rotI, indexConfTreeEmat_J, rotJ, currentE - toAddSubtract);

                        }
                    }
                }
            }
        }
    }

    private void updateSubsetPruneMat(PruningMatrix currentPruneMat, PruningMatrix seqNodePruneMat, ArrayList<Integer> subsetOfPositions) {
        Collections.sort(subsetOfPositions);
        int numPos = subsetOfPositions.size();
        for (int i = 0; i < numPos; i++) {
            int origPosNum = subsetOfPositions.get(i);
            for (int rotI = 0; rotI < currentPruneMat.oneBody.get(i).size(); rotI++) {
                currentPruneMat.setOneBody(i, rotI, seqNodePruneMat.getOneBody(origPosNum, rotI));
            }
            for (int j = i + 1; j < numPos; j++) {
                int origPosNum_J = subsetOfPositions.get(j);
                int origPosNum_I = subsetOfPositions.get(i);
                for (int rotJ = 0; rotJ < currentPruneMat.pairwise.get(j).get(i).size(); rotJ++) {
                    for (int rotI = 0; rotI < currentPruneMat.pairwise.get(j).get(i).get(rotJ).size(); rotI++) {
                        currentPruneMat.setPairwise(j, rotJ, i, rotI, seqNodePruneMat.getPairwise(origPosNum_J, rotJ, origPosNum_I, rotI));
                    }
                }
            }
        }
    }

    private void updateEnergyMatrix_CrossTerms(SearchProblemSuper searchProblem, ArrayList<Integer> subsetOfPositions, boolean[][] interactionGraph, boolean negateEnergies) {
        EnergyMatrix origEmat = searchProblem.emat;
        PruningMatrix origPruneMat = searchProblem.pruneMat;

        EnergyMatrix updatedEmat = new EnergyMatrix(searchProblem.confSpaceSuper, searchProblem.emat.getPruningInterval());

        int numPos = searchProblem.confSpaceSuper.numPos;
        int[] numRotPerPos = searchProblem.confSpaceSuper.getNumRCsAtPos();

        // Now we can update the energy matrix
        // For cross-terms we only add pairwise energies
        for (int i = 0; i < numPos; i++) {
            int origPosNum_I = subsetOfPositions.get(i);
            //Set all one-body terms to 0.0
            for (int rotI = 0; rotI < numRotPerPos[i]; rotI++) {
                updatedEmat.setOneBody(i, rotI, 0.0);
            }
            for (int j = i + 1; j < numPos; j++) {
                int origPosNum_J = subsetOfPositions.get(j);
                // check if they are interacting
                if (interactionGraph[i][j]) {
                    // they interact, so add their pairwise terms
                    for (int rotI = 0; rotI < numRotPerPos[i]; rotI++) {
                        for (int rotJ = 0; rotJ < numRotPerPos[j]; rotJ++) {
                            double pairwiseE = origEmat.getPairwise(origPosNum_J, rotJ, origPosNum_I, rotI);
                            if (negateEnergies) {
                                pairwiseE = -pairwiseE;
                            }
                            boolean pairIsPruned = origPruneMat.getPairwise(origPosNum_I, rotI, origPosNum_J, rotJ);
                            if (pairIsPruned) {
                                updatedEmat.setPairwise(j, rotJ, i, rotI, Double.POSITIVE_INFINITY);
                            } else {
                                updatedEmat.setPairwise(j, rotJ, i, rotI, pairwiseE);
                            }
                        }
                    }
                } else {
                    for (int rotI = 0; rotI < numRotPerPos[i]; rotI++) {
                        for (int rotJ = 0; rotJ < numRotPerPos[j]; rotJ++) {
                            updatedEmat.setPairwise(j, rotJ, i, rotI, 0.0);
                        }
                    }
                }
            }
        }
        searchProblem.emat = updatedEmat;
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
                int newPosNum_J = allPositions.indexOf(posJ);
                interactionGraph[newPosNum_I][newPosNum_J] = true;
                interactionGraph[newPosNum_J][newPosNum_I] = true;
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
    private ArrayList<Integer> getLigandAssignedPosNums(COMETSNodeSuper seqNode, boolean useBoundState) {
        int[] partialSeq = seqNode.getNodeAssignments();
        //Get the mapping between mutable position in sequence node assignment and
        //the corresponding position number in the confSpace for bound or unbound ligand
        ArrayList<Integer> mutable2PosNum = useBoundState ? this.mutable2StatePosNums.get(0) : this.mutable2StatePosNums.get(1);
        //Get the corresponding searchProblem for bound or unbound ligand
        SearchProblemSuper searchProblem = useBoundState ? this.mutableSearchProblems[0] : this.mutableSearchProblems[1];

        ArrayList<Integer> ligandAssignedPosNums = new ArrayList<>();

        //Iterate over flexible position numbers 
        for (int posNum = 0; posNum < searchProblem.confSpaceSuper.numPos; posNum++) {

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
    private ArrayList<Integer> getLigandUnassignedPosNums(COMETSNodeSuper seqNode, boolean useBoundState) {
        int[] partialSeq = seqNode.getNodeAssignments();

        //Get the mapping between mutable position in sequence node assignment and
        //the corresponding position number in the confSpace for bound or unbound ligand
        ArrayList<Integer> mutable2PosNum = useBoundState ? this.mutable2StatePosNums.get(0) : this.mutable2StatePosNums.get(1);
        //Get the corresponding searchProblem for bound or unbound ligand
        SearchProblemSuper searchProblem = useBoundState ? this.mutableSearchProblems[0] : this.mutableSearchProblems[1];

        ArrayList<Integer> ligandUnassignedPosNums = new ArrayList<>();

        //Iterate over flexible position numbers
        for (int posNum = 0; posNum < searchProblem.confSpaceSuper.numPos; posNum++) {

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
        SearchProblemSuper searchSpace = useBoundState ? mutableSearchProblems[0] : nonMutableSearchProblem;

        ArrayList<Integer> proteinPosNums = new ArrayList<>();

        for (int posNum = 0; posNum < searchSpace.confSpaceSuper.numPos; posNum++) {
            //If we are using the bound state we must check if posNum
            //belongs to protein or ligand
            if (useBoundState) {
                //Check if pos belongs to nonmutable (protein) strand
                if (!this.boundResNumToIsMutableStrand.get(posNum)) {
                    proteinPosNums.add(posNum);
                }
            } else {//If we are using the unbound state, every posNum is part of protein
                proteinPosNums.add(posNum);
            }
        }
        return proteinPosNums;
    }
    
        private String[] getSequence(COMETSNodeSuper node){
        int[] assignments = node.getNodeAssignments();
        int numMotPos = assignments.length;
        String[] sequence = new String[numMotPos];

        for (int mutPos = 0; mutPos<numMotPos; mutPos++){
            int aaTypeVal = assignments[mutPos];
            String aaType;
            if (aaTypeVal==-1){
                aaType = "XXX";
            }
            else{
                aaType = this.AATypeOptions.get(mutPos).get(aaTypeVal);
            }
            sequence[mutPos] = aaType;
        }
        return sequence;
    }
    
    private void printSequence(String[] sequence){
        StringBuffer buffer = new StringBuffer();
        for (String aaType : sequence){
            buffer.append(" "+aaType);
        }
        System.out.println(buffer);
    }
}
