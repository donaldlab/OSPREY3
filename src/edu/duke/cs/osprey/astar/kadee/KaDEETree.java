/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.kadee;

import edu.duke.cs.osprey.astar.comets.*;
import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.AStarTree;
import edu.duke.cs.osprey.astar.ConfTreeSuper;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.SuperRCTuple;
import edu.duke.cs.osprey.confspace.SearchProblemSuper;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PrunerSuper;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.stream.Collectors;
import org.apache.commons.lang.ArrayUtils;

/**
 * This class behaves similarly to COMETSTree and COMETSTreeSuper It optimizes
 * for the sequences with the best K* score
 *
 * @author hmn5
 */
public class KaDEETree extends AStarTree {

    int numTreeLevels; //number of mutable residues

    LME objFcn; //objective function to minimize
    LME[] constraints; //constraints on our sequence

    ArrayList<ArrayList<String>> AATypeOptions; //The allowed amino-acids at each level

    int numMaxMut; //number of mutatations allowed away from wtSeq (-1 means no cap)
    String wtSeq[]; //wt sequence
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

    //Maps the bound res num to the corresponding unbound emat
    HashMap<Integer, EnergyMatrix> boundResNumToUnboundEmat;
    //Maps the bound res num to the corresponding unbound res num
    HashMap<Integer, Integer> boundResNumToUnboundResNum;
    //Maps the bound res num to boolean that is true if res num is part of mutable
    //strand
    HashMap<Integer, Boolean> boundResNumToIsMutableStrand;
    //determines if two residues are on the same strand
    boolean[][] belongToSameStrand;

    public KaDEETree(int numTreeLevels, LME objFcn, LME[] constraints,
            ArrayList<ArrayList<String>> AATypeOptions, int numMaxMut, String[] wtSeq,
            int numStates, SearchProblemSuper[] stateSP,
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
        KaDEENode seqNode = (KaDEENode) curNode;
        ArrayList<AStarNode> ans = new ArrayList<>();

        if (seqNode.isFullyDefined()) {
            ans.add(seqNode);
            return ans;
        } else {
            //expand next position...
            int[] curAssignments = seqNode.getNodeAssignments();

            for (int splitPos = 0; splitPos < numTreeLevels; splitPos++) {
                if (curAssignments[splitPos] < 0) {//we can split this level

                    for (int aa = 0; aa < AATypeOptions.get(splitPos).size(); aa++) {

                        int[] childAssignments = curAssignments.clone();
                        childAssignments[splitPos] = aa;
                        
                        UpdatedPruningMatrixSuper[] childPruneMat = new UpdatedPruningMatrixSuper[numStates];
                        for (int state = 0; state<numStates; state++){
                            childPruneMat[state] = doChildPruning(state, seqNode.pruneMat[state], splitPos, aa);
                        }
                        
                        KaDEENode childNode = new KaDEENode(childAssignments, childPruneMat);
                        
                        if (splitPos == numTreeLevels-1){//sequence is fully defined...make conf trees
                            makeSeqConfTrees(childNode);
                        }

                        //TODO: create boundFreeEnergyChange()
                        //childNode.setScore(boundFreeEnergyChange(childNode));
                        
                        ans.add(childNode);
                    }
                    
                    return ans;
                }
            }
            
            throw new RuntimeException("ERROR: Not splittable position found but sequence not fully defined...");
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
        return ans;
    }

     private void makeSeqConfTrees(KaDEENode node) {
        //Given a node with a fully defined sequence, build its conformational search trees
        //for each state
        //If a state has no viable conformations, leave it null, with stateUB[state] = inf

        node.stateTrees = new ConfTreeSuper[numStates];

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

                node.stateTrees[state].initQueue(rootNode);//allocate queue and add root node
            } else {//no confs available for this state!
                node.stateTrees[state] = null;
            }
        }
    }   
    
    
    
    //For flexible position in bound matrix (0,1,2,...,numRes-1) we map to the corresponding
    //unbound energy matrix
    public HashMap<Integer, EnergyMatrix> getBoundPosNumToUnboundEmat() {
        SearchProblemSuper boundState = this.mutableSearchProblems[0];
        //Get res number from each flexible position in the bound state
        List<Integer> resNumsBound = boundState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound mutable state
        SearchProblemSuper unBoundMutableState = this.mutableSearchProblems[1];
        List<Integer> resNumsUnboundMutable = unBoundMutableState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound non-mutable state
        SearchProblemSuper unBoundNonMutableState = this.nonMutableSearchProblem;
        List<Integer> resNumsUnboundNonMutable = unBoundNonMutableState.confSpaceSuper.posFlexSuper.stream()
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
        List<Integer> resNumsBound = boundState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound mutable state
        SearchProblemSuper unBoundMutableState = this.mutableSearchProblems[1];
        List<Integer> resNumsUnboundMutable = unBoundMutableState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound non-mutable state
        SearchProblemSuper unBoundNonMutableState = this.nonMutableSearchProblem;
        List<Integer> resNumsUnboundNonMutable = unBoundNonMutableState.confSpaceSuper.posFlexSuper.stream()
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
        List<Integer> resNumsBound = boundState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toList());

        //Get res number for each flexible position in the unbound mutable state
        SearchProblemSuper unBoundMutableState = this.mutableSearchProblems[1];
        List<Integer> resNumsUnboundMutable = unBoundMutableState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toList());

        //Get res number for each flexible position in the unbound non-mutable state
        SearchProblemSuper unBoundNonMutableState = this.nonMutableSearchProblem;
        List<Integer> resNumsUnboundNonMutable = unBoundNonMutableState.confSpaceSuper.posFlexSuper.stream()
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
        List<Integer> resNumsBound = boundState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound mutable state
        SearchProblemSuper unBoundMutableState = this.mutableSearchProblems[1];
        List<Integer> resNumsUnboundMutable = unBoundMutableState.confSpaceSuper.posFlexSuper.stream()
                .map(posFlex -> posFlex.resList.get(0).resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        //Get res number for each flexible position in the unbound non-mutable state
        SearchProblemSuper unBoundNonMutableState = this.nonMutableSearchProblem;
        List<Integer> resNumsUnboundNonMutable = unBoundNonMutableState.confSpaceSuper.posFlexSuper.stream()
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
}
