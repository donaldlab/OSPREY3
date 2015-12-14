/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.iminmsd;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.AStarTree;
import edu.duke.cs.osprey.astar.ConfTreeSuper;
import edu.duke.cs.osprey.astar.comets.LME;
import edu.duke.cs.osprey.astar.comets.UpdatedPruningMatrixSuper;
import edu.duke.cs.osprey.confspace.PositionConfSpaceSuper;
import edu.duke.cs.osprey.confspace.SearchProblemSuper;
import edu.duke.cs.osprey.confspace.SuperRCTuple;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.pruning.Pruner;
import edu.duke.cs.osprey.pruning.PrunerSuper;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;
import org.apache.commons.lang.ArrayUtils;

/**
 * This class is an extension of COMETS for continuous Multi-state Design (MSD)
 * within the iMinDEE framework for Bound vs Unbound MSD
 *
 * @author hmn5
 */
public class iMinMSDTree extends AStarTree {

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

    SearchProblemSuper nonMutableSearchProblem;

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
    HashMap<Integer, Boolean> boundResNumToIsLigand;
    //determines if two residues are on the same strand
    boolean[][] belongToSameStrand;

    EnergyMatrix maxInterfaceBoundEmat;

    public iMinMSDTree(int numTreeLevels, LME objFcn, LME[] constraints,
            ArrayList<ArrayList<String>> AATypeOptions, int numMaxMut, String[] wtSeq,
            int numStates, SearchProblemSuper[] stateSP, SearchProblemSuper nonMutableState,
            ArrayList<ArrayList<Integer>> mutable2StatePosNums) {

        this.numTreeLevels = numTreeLevels;
        this.objFcn = objFcn;
        this.constraints = constraints;
        this.AATypeOptions = AATypeOptions;
        this.numMaxMut = numMaxMut;
        this.wtSeq = wtSeq;
        this.numStates = numStates;
        this.mutableSearchProblems = stateSP;
        this.nonMutableSearchProblem = nonMutableState;
        this.mutable2StatePosNums = mutable2StatePosNums;

        stateNumPos = new int[numStates];
        for (int state = 0; state < numStates; state++) {
            stateNumPos[state] = stateSP[state].confSpaceSuper.numPos;
        }

        this.boundResNumToUnboundEmat = getBoundPosNumToUnboundEmat();
        this.boundResNumToUnboundResNum = getBoundPosNumToUnboundPosNum();
        this.boundResNumToIsLigand = getBoundPosNumberToIsMutableStrand();
        this.belongToSameStrand = getSameStrandMatrix();

        this.maxInterfaceBoundEmat = calcMaxInterfaceEmat();

        //this.objFcn.setConstTerm(-calcGMEC(nonMutableState));
        //TODO: Fix this 
        this.setNumCores(8);
    }

    /*
    @Override
    public ArrayList<AStarNode> getChildren(AStarNode curNode) {
        iMinMSDNode seqNode = (iMinMSDNode) curNode;
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
                        for (int state = 0; state < numStates; state++) {
                            childPruneMat[state] = doChildPruning(state, seqNode.pruneMat[state], splitPos, aa);
                        }

                        iMinMSDNode childNode = new iMinMSDNode(childAssignments, childPruneMat);

                        if (splitPos == numTreeLevels - 1) {//sequence is fully defined...make conf trees
                            makeSeqConfTrees(childNode);
                        }

                        childNode.setScore(boundLME(childNode));
                        ans.add(childNode);
                    }

                    return ans;
                }
            }

            throw new RuntimeException("ERROR: Not splittable position found but sequence not fully defined...");
        }
    }
*/
    @Override
    public ArrayList<AStarNode> getChildren(AStarNode curNode) {
        iMinMSDNode seqNode = (iMinMSDNode) curNode;
        ArrayList<AStarNode> ans = new ArrayList<>();

        if (seqNode.isFullyDefined()) {
            ans.add(seqNode);
            return ans;
        } else {
            //expand next position...

            ArrayList<iMinMSDNode> children = new ArrayList<>();
            int[] curAssignments = seqNode.getNodeAssignments();

            for (int splitPos = 0; splitPos < numTreeLevels; splitPos++) {
                if (curAssignments[splitPos] < 0) {//we can split this level

                    ArrayList<Integer> aaTypes = new ArrayList<>();
                    for (int aa = 0; aa < AATypeOptions.get(splitPos).size(); aa++) {

                        int[] childAssignments = curAssignments.clone();
                        childAssignments[splitPos] = aa;

                        UpdatedPruningMatrixSuper[] childPruneMat = new UpdatedPruningMatrixSuper[numStates];
                        for (int state = 0; state < numStates; state++) {
                            childPruneMat[state] = this.doChildPruning(state, seqNode.pruneMat[state], splitPos, aa);
                        }

                        iMinMSDNode childNode = new iMinMSDNode(childAssignments, childPruneMat);

                        if (splitPos == numTreeLevels - 1) {//sequence is fully defined...make conf trees
                            makeSeqConfTrees(childNode);
                        }
                        children.add(childNode);
                    }
                    children.parallelStream().forEach(node -> {
                        try {
                            node.setScore(boundLME(node));
                            ans.add(node);
                        } catch (Exception e) {
                            System.out.println(e.getMessage());
                            System.exit(1);
                        }
                    });
                    return ans;
                }
            }

            throw new RuntimeException("ERROR: Not splittable position found but sequence not fully defined...");
        }
    }

    private double boundLME(iMinMSDNode seqNode) {
        if (seqNode.isFullyDefined()) {
            double score = 0.0;
            //TODO: Calc MAP for bound and unbound
            //Add objFcn.constTerm
            return score;
        } else {
            double score = maxInterfaceBound(seqNode);
            //TODO: Max interface bound + constTerms
            return score;
        }
    }

    private double maxInterfaceBound(iMinMSDNode seqNode) {
        SearchProblemSuper boundSP = mutableSearchProblems[0];
        SearchProblemSuper ligandSP = mutableSearchProblems[1];

        ArrayList<Integer> proteinPosNums = getProteinPosNums(true);
        ArrayList<Integer> ligandPosNums = getLigandAssignedPosNums(seqNode, true);
        ligandPosNums.addAll(getLigandUnassignedPosNums(seqNode, true));
        Collections.sort(ligandPosNums);

        ArrayList<Integer> allPos = new ArrayList<>();
        allPos.addAll(ligandPosNums);
        allPos.addAll(proteinPosNums);
        Collections.sort(allPos);

        //SearchProblemSuper sp_Interface = boundSP.getSubsetSearchProblem(allPos);
        SearchProblemSuper sp_Interface = (SearchProblemSuper) ObjectIO.deepCopy(boundSP);
        sp_Interface.emat.setConstTerm(0.0);
        updateSubsetPruneMat(sp_Interface.pruneMat, seqNode.pruneMat[0], allPos);
        subtractInternalEnergies(sp_Interface.emat, boundSP.emat, allPos, ligandPosNums);
        subtractPairwiseEnergies(sp_Interface.emat, boundSP.emat, allPos, ligandPosNums);
        setMaxInterfaceLigandIntraEnergy(sp_Interface.emat);
        
        sp_Interface.fullConfE = genMaxInterfaceEnergyFunction(sp_Interface);
               
        double lowerBound = calcGMEC(sp_Interface);
        /*
        ConfTreeSuper confTree = new ConfTreeSuper(sp_Interface);
        double lowerBound = confTree.energyNextConf();
        */

        return lowerBound + objFcn.getConstTerm();
    }

    private double fullSequenceBound(iMinMSDNode seqNode) {
        SearchProblemSuper boundSP = mutableSearchProblems[0];
        SearchProblemSuper ligandSP = mutableSearchProblems[1];

        ArrayList<Integer> proteinPosNums = getProteinPosNums(true);
        ArrayList<Integer> ligandPosNums = getLigandAssignedPosNums(seqNode, true);
        ligandPosNums.addAll(getLigandUnassignedPosNums(seqNode, true));
        Collections.sort(ligandPosNums);

        ArrayList<Integer> allPos = new ArrayList<>();
        allPos.addAll(ligandPosNums);
        allPos.addAll(proteinPosNums);
        Collections.sort(allPos);

        SearchProblemSuper spBound = boundSP.getSubsetSearchProblem(allPos);
        updateSubsetPruneMat(spBound.pruneMat, seqNode.pruneMat[0], allPos);
        double gmecBound = calcGMEC(spBound);

        ArrayList<Integer> ligandUnbound = new ArrayList<>();
        ligandUnbound.addAll(getLigandAssignedPosNums(seqNode, false));
        ligandUnbound.addAll(getLigandUnassignedPosNums(seqNode, false));
        Collections.sort(ligandUnbound);

        SearchProblemSuper spUnbound = ligandSP.getSubsetSearchProblem(ligandUnbound);
        updateSubsetPruneMat(spUnbound.pruneMat, seqNode.pruneMat[1], ligandUnbound);
        double gmecUnbound = calcGMEC(spUnbound);

        double score = gmecBound - gmecUnbound;
        return 0.0;
    }

    /**
     * Recalculates the intra-energy terms for the searchProblem so that the
     * ligand intra-shell term contains just ligand intra and protein template
     * terms
     *
     * @param sp the searchProblem
     */
    private void recalculateLigandIntraTerms(SearchProblemSuper sp) {
        ArrayList<Residue> boundShellResProtein = getBoundShellResiduesIfProtein(sp);
        EnergyMatrixCalculator ematCalc = new EnergyMatrixCalculator(sp.confSpaceSuper, boundShellResProtein);
        ematCalc.reCalculateIntraTerms(sp.emat, sp.pruneMat, boundResNumToIsLigand);
    }

    /*
     private EnergyFunction getMaxInterfaceEnergyFunction(SearchProblemSuper maxInterSearchProb, HashMap<Integer, Boolean> boundPosNum2IsLigandStrand){
     SearchProblemSuper ligandSP = mutableSearchProblems[1];
     SearchProblemSuper proteinSP = nonMutableSearchProblem;
     }
     */
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

    private void makeSeqConfTrees(iMinMSDNode node) {
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

    @Override
    public boolean isFullyAssigned(AStarNode node) {
        //HMN: TODO: We need to decide when a K* score is fully calculated

        if (!node.isFullyDefined()) {
            return false;
        }
        return true;
    }

    @Override
    public AStarNode rootNode() {
        int[] conf = new int[numTreeLevels];
        Arrays.fill(conf, -1);//indicates the sequence is not assigned

        PruningMatrix[] pruneMats = new PruningMatrix[numStates];
        for (int state = 0; state < numStates; state++) {
            pruneMats[state] = mutableSearchProblems[state].pruneMat;
        }

        iMinMSDNode root = new iMinMSDNode(conf, pruneMats);
        //TODO: root.setScore(boundLME(root,objFcn));
        return root;
    }

    @Override
    public boolean canPruneNode(AStarNode node) {
        //Check if node can be pruned
        //This is traditionally based on constraints, thought we could pruned nodes
        //that are provably suboptimal

        //TODO: Implement constraints as in COMETS if desired
        iMinMSDNode seqNode = (iMinMSDNode) node;

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
        if (false) {
            for (LME constr : constraints) {
                if ((constr.getCoeffs()[0] == 1.0) && (constr.getCoeffs()[1] == -1.0)) {
                    if (maxInterfaceBound(seqNode) + constr.getConstTerm() > 0) {
                        return true;
                    }
                } else if ((constr.getCoeffs()[0] == 1.0) && (constr.getCoeffs()[1] == 0.0)) {
                    UpdatedPruningMatrixSuper pruneMatBound = new UpdatedPruningMatrixSuper(seqNode.pruneMat[0]);
                    PrunerSuper dee = new PrunerSuper(mutableSearchProblems[0], pruneMatBound, false, 100, 0.0, false, false, false);
                    int oldNumUpdates = 0;
                    int numUpdates = 0;
                    do {
                        oldNumUpdates = numUpdates;
                        dee.prune("GOLDSTEIN");
                        numUpdates = pruneMatBound.countUpdates();
                    } while (numUpdates > oldNumUpdates);

                    ConfTreeSuper confTreeBound = new ConfTreeSuper(mutableSearchProblems[0], pruneMatBound, false);
                    double gmecEBound = confTreeBound.energyNextConf();
                    if (gmecEBound + constr.getConstTerm() > 0) {
                        return true;
                    }
                } else if ((constr.getCoeffs()[0] == 0.0) && (constr.getCoeffs()[1] == 1.0)) {
                    UpdatedPruningMatrixSuper pruneMatUnbound = new UpdatedPruningMatrixSuper(seqNode.pruneMat[1]);
                    PrunerSuper dee = new PrunerSuper(mutableSearchProblems[1], pruneMatUnbound, false, 100, 0.0, false, false, false);
                    int oldNumUpdates = 0;
                    int numUpdates = 0;
                    do {
                        oldNumUpdates = numUpdates;
                        dee.prune("GOLDSTEIN");
                        numUpdates = pruneMatUnbound.countUpdates();
                    } while (numUpdates > oldNumUpdates);

                    ConfTreeSuper confTreeUnbound = new ConfTreeSuper(mutableSearchProblems[1], pruneMatUnbound, false);
                    double gmecEUnbound = confTreeUnbound.energyNextConf();
                    if (gmecEUnbound + constr.getConstTerm() > 0) {
                        return true;
                    }
                }
            }
        }
        return false;
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
    private ArrayList<Integer> getLigandAssignedPosNums(iMinMSDNode seqNode, boolean useBoundState) {
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
            if (this.boundResNumToIsLigand.get(posNum) || !(useBoundState)) {
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
    private ArrayList<Integer> getLigandUnassignedPosNums(iMinMSDNode seqNode, boolean useBoundState) {
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
            if (this.boundResNumToIsLigand.get(posNum) || !(useBoundState)) {
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
                if (!this.boundResNumToIsLigand.get(posNum)) {
                    proteinPosNums.add(posNum);
                }
            } else {//If we are using the unbound state, every posNum is part of protein
                proteinPosNums.add(posNum);
            }
        }
        return proteinPosNums;
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

    public ArrayList<Residue> getBoundShellResiduesIfProtein(SearchProblemSuper sp) {
        SearchProblemSuper proteinSP = nonMutableSearchProblem;
        ArrayList<Residue> boundShellResIsProtein = new ArrayList<>();
        for (Residue shellRes : sp.shellResidues) {
            for (Residue proteinShellRes : proteinSP.shellResidues) {
                if (shellRes.resNum == proteinShellRes.resNum) {
                    boundShellResIsProtein.add(shellRes);
                    break;
                }
            }
        }
        return boundShellResIsProtein;
    }

    private EnergyFunction genMaxInterfaceEnergyFunction(SearchProblemSuper boundSP) {
        ArrayList<Residue> boundShellResProtein = getBoundShellResiduesIfProtein(boundSP);

        //Separate flexible residues into those that belong to protein and ligand
        ArrayList<Residue> flexibleLigandResidues = new ArrayList<>();
        ArrayList<Residue> flexibleProteinResidues = new ArrayList<>();

        for (int pos = 0; pos < boundSP.confSpaceSuper.numPos; pos++) {
            PositionConfSpaceSuper pcs = boundSP.confSpaceSuper.posFlexSuper.get(pos);
            if (pcs.resList.size() != 1) {
                throw new RuntimeException("ERROR: Max Interface Energy Bound currently only supports one res per pos");
            }
            Residue flexRes = pcs.resList.get(0);
            if (boundResNumToIsLigand.get(pos)) {
                flexibleLigandResidues.add(flexRes);
            } else {
                flexibleProteinResidues.add(flexRes);
            }
        }
        EnergyFunctionGenerator efg = EnvironmentVars.curEFcnGenerator;
        MultiTermEnergyFunction maxInterfaceEFunc = new MultiTermEnergyFunction();

        for (Residue ligandRes : flexibleLigandResidues) {
            for (Residue proteinRes : flexibleProteinResidues) {
                EnergyFunction pairE = efg.resPairEnergy(ligandRes, proteinRes);
                maxInterfaceEFunc.addTerm(pairE);
            }
        }

        for (int i = 0; i < flexibleProteinResidues.size(); i++) {
            for (int j = i + 1; j < flexibleProteinResidues.size(); j++) {
                Residue proteinResI = flexibleProteinResidues.get(i);
                Residue proteinResJ = flexibleProteinResidues.get(j);
                EnergyFunction pairE = efg.resPairEnergy(proteinResI, proteinResJ);
                maxInterfaceEFunc.addTerm(pairE);
            }
        }
        //ligand protein residues only interact with protein shell residues for 
        //max interface bound
        for (Residue ligandRes : flexibleLigandResidues) {
            EnergyFunction oneBodyE = efg.singleResEnergy(ligandRes);
            maxInterfaceEFunc.addTerm(oneBodyE);
            for (Residue proteinShellRes : boundShellResProtein) {
                EnergyFunction twoBodyE = efg.resPairEnergy(ligandRes, proteinShellRes);
                maxInterfaceEFunc.addTerm(twoBodyE);
            }
        }

        //protein flexible residues interact with all shell residues
        for (Residue proteinRes : flexibleProteinResidues) {
            EnergyFunction oneBodyE = efg.singleResEnergy(proteinRes);
            maxInterfaceEFunc.addTerm(oneBodyE);
            for (Residue shellRes : boundSP.shellResidues) {
                EnergyFunction twoBodyE = efg.resPairEnergy(proteinRes, shellRes);
                maxInterfaceEFunc.addTerm(twoBodyE);
            }
        }

        return maxInterfaceEFunc;
    }

    private double calcGMEC(SearchProblemSuper searchProb) {

        double IO = Double.POSITIVE_INFINITY;
        double Ew = 0;
        double lowestBound = Double.POSITIVE_INFINITY;

        int[] GMECConf = null;
        double bestESoFar = Double.POSITIVE_INFINITY;
        boolean needToRepeat;

        do {
            needToRepeat = false;

            int numUpdates = 0;
            int oldNumUpdates;
            UpdatedPruningMatrixSuper pruneMat = new UpdatedPruningMatrixSuper(searchProb.pruneMat);
            PrunerSuper dee = new PrunerSuper(searchProb, pruneMat, false, 100, IO, false, false, false);
            do {
                oldNumUpdates = numUpdates;
                dee.prune("GOLDSTEIN");
                numUpdates = pruneMat.countUpdates();
            } while (numUpdates > oldNumUpdates);

            ConfTreeSuper confTree = new ConfTreeSuper(searchProb, pruneMat, false);

            double lowerBound;
            int conformationCount = 0;

            System.out.println();
            System.out.println("BEGINNING CONFORMATION ENUMERATION");
            System.out.println();

            long confSearchStartTime = System.currentTimeMillis();

            do {
                int conf[] = confTree.nextConf();

                if (conf == null) {
                    lowerBound = Double.POSITIVE_INFINITY;
                } else {
                    double confE = searchProb.minimizedEnergy(conf);
                    //System.out.println("Actual Energy: "+confE);
                    lowerBound = searchProb.lowerBound(conf);

                    System.out.println("Lower Bound: " + lowerBound + " Minimized Energy: " + confE);
                    System.out.println();

                    if (confE < bestESoFar) {
                        bestESoFar = confE;
                        GMECConf = conf;
                    }

                    lowestBound = Math.min(lowestBound, lowerBound);
                }

                if ((lowerBound > lowestBound + IO) && (bestESoFar > lowestBound + IO)) {
                    IO = bestESoFar - lowerBound + 0.001;
                    needToRepeat = true;
                    break;
                }
                if (bestESoFar == Double.POSITIVE_INFINITY) {
                    System.out.println("A* returned no conformations.");
                    break;
                }
            } while (bestESoFar >= lowerBound);
            double confSearchTimeMinutes = (System.currentTimeMillis() - confSearchStartTime) / 60000.0;
            System.out.println("Conf search time (minutes): " + confSearchTimeMinutes);
        } while (needToRepeat);
        System.out.println("GMEC energy: " + bestESoFar);
        return bestESoFar;

    }

    private EnergyMatrix calcMaxInterfaceEmat() {
        SearchProblemSuper boundSP = (SearchProblemSuper) ObjectIO.deepCopy(this.mutableSearchProblems[0]);
        recalculateLigandIntraTerms(boundSP);
        return boundSP.emat;
    }

    private void setMaxInterfaceLigandIntraEnergy(EnergyMatrix emat) {
        int numPos = emat.oneBody.size();
        emat.oneBody = new ArrayList<>();
        for (int pos = 0; pos < numPos; pos++) {
            emat.oneBody.add(this.maxInterfaceBoundEmat.oneBody.get(pos));
        }
    }

    public void setNumCores(int n) {
        int cores = Math.min(n, Runtime.getRuntime().availableProcessors());
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism",
                String.valueOf(cores));
    }
}
