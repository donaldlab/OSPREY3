//package edu.duke.cs.osprey.astar.ewakstar;
//
//import edu.duke.cs.osprey.astar.AStarTree;
//import edu.duke.cs.osprey.astar.ConfTree;
//import edu.duke.cs.osprey.astar.FullAStarNode;
//import edu.duke.cs.osprey.astar.comets.UpdatedPruningMatrix;
//import edu.duke.cs.osprey.confspace.*;
//import edu.duke.cs.osprey.ematrix.EnergyMatrix;
//import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
//import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
//import edu.duke.cs.osprey.gmec.PrecomputedMatrices;
//import edu.duke.cs.osprey.pruning.PruningMatrix;
//
//import java.util.ArrayList;
//import java.util.Arrays;
//
//public class NewEWAKStarTree extends AStarTree<EWAKStarNode> {
//
//    private int numTreeLevels;//number of residues with sequence changes
//
//    private int numMutable;
//    private boolean useExact;
//    private ArrayList<ArrayList<String>> AATypeOptions;
//    //COMETreeNode.assignments assigns each level an index in AATypeOptions.get(level), and thus an AA type
//    //If -1, then no assignment yet
//
//    private Sequence wtSeq;
//    private String[] wtSeqList;
//    //information on states
//
//    //description of each state
//    private SimpleConfSpace confSpace;//their conformational spaces
//    private EnergyMatrix emat;
//    private PruningMatrix pmat;
//
//
//    private ArrayList<Integer> mutablePosNums;
//    //mutable2StatePosNum.get(state) maps levels in this tree to flexible positions for state
//    //(not necessarily an onto mapping)
//
//    private int stateNumPos;
//    double wtMinimizedEnergy = Double.POSITIVE_INFINITY;
//
//
//    ConfEnergyCalculator confECalc = null;//only needed if we want minimized structs.  one per state like the other arrays
//
//    public NewEWAKStarTree(boolean useExact, int numMutable, int numTreeLevels, ArrayList<ArrayList<String>> AATypeOptions,
//                           Sequence wtSeq, SimpleConfSpace confSpace,
//                           EnergyMatrix emat, PruningMatrix pmat, ArrayList<Integer> mutablePosNums,
//                         ConfEnergyCalculator confECalc) {
//
//        this.useExact = useExact;
//        this.numMutable = numMutable;
//        this.numTreeLevels = numTreeLevels;
//        this.AATypeOptions = AATypeOptions;
//        this.wtSeq = wtSeq;
//        this.wtSeqList = NewEWAKStarDoer.toEWAKStarString(wtSeq).split("_");
//        this.confSpace = confSpace;
//        this.emat = emat;
//        this.pmat = pmat;
//        this.mutablePosNums = mutablePosNums;
//        this.confECalc = confECalc;
//
//        stateNumPos = confSpace.getNumPos();
//    }
//
//    public void setAATypeOptions(ArrayList<ArrayList<String>> newAATypeOptions){
//        this.AATypeOptions = newAATypeOptions;
//    }
//
//    @Override
//    public ArrayList<EWAKStarNode> getChildren(EWAKStarNode seqNode) {
//
//        ArrayList<EWAKStarNode> seqNodeChildren = new ArrayList<>();
//
//        if(seqNode.isFullyDefined()){
//            seqNode.expandConfTree();
//            seqNode.setScore( boundLME(seqNode) );
//            seqNodeChildren.add(seqNode);
//            return seqNodeChildren;
//        }
//        else{
//            //expand next position...
//            int curAssignments[] = seqNode.getNodeAssignments();
//
//            for(int splitPos=0; splitPos<numTreeLevels; splitPos++){
//                if(curAssignments[splitPos] < 0){//can split this level
//
//                    for(int aa=0; aa<AATypeOptions.get(splitPos).size(); aa++){
//
//                        int childAssignments[] = curAssignments.clone();
//                        childAssignments[splitPos] = aa;
//
//                        UpdatedPruningMatrix childPruneMat = doChildPruning(seqNode.pruneMat, splitPos, aa);
//
//                        EWAKStarNode childNode = new EWAKStarNode(childAssignments, childPruneMat);
//
//                        if(splitPos==numTreeLevels-1){//sequence now fully defined...make conf trees
//                            makeSeqConfTrees(childNode);
//                        }
//
//                        childNode.setScore( boundLME(childNode) );
//                        seqNodeChildren.add(childNode);
//                    }
//
//                    return seqNodeChildren;
//                }
//            }
//
//            throw new RuntimeException("ERROR: No splittable position found but sequence not fully defined...");
//        }
//
//    }
//
//
//    private UpdatedPruningMatrix doChildPruning(PruningMatrix parentMat, int splitPos, int aa){
//        //Create an update to parentMat (without changing parentMat)
//        //to reflect that splitPos has been assigned an amino-acid type
//
//        String assignedAAType = AATypeOptions.get(splitPos).get(aa);
//
//        UpdatedPruningMatrix ans = new UpdatedPruningMatrix(parentMat);
//        int posAtState = mutablePosNums.get(splitPos);
//
//        //first, prune all other AA types at splitPos
//        for(int rc : parentMat.unprunedRCsAtPos(posAtState)){
//            String rcAAType = confSpace.positions.get(posAtState).resConfs.get(rc).template.name;
//
//            if( ! rcAAType.equalsIgnoreCase(assignedAAType)){
//                ans.markAsPruned(new RCTuple(posAtState,rc));
//            }
//        }
//
////        removed pruning done in COMETs - not necessary and just takes time here.
//
//        return ans;
//    }
//
//
//
//    private void makeSeqConfTrees(EWAKStarNode node){
//        //Given a node with a fully defined sequence, build its conformational search trees
//        //If a state has no viable conformations, leave it null, with stateUB[state] = inf
//
//        //first make sure there are RCs available at each position
//        boolean RCsAvailable = true;
//        for(int pos=0; pos<stateNumPos; pos++){
//            if(node.pruneMat.unprunedRCsAtPos(pos).isEmpty()){
//                RCsAvailable = false;
//                break;
//            }
//        }
//
//        if(RCsAvailable) {
//            node.stateTree = new ConfTree<>(
//                    new FullAStarNode.Factory(stateNumPos),
//                    false,
//                    emat,
//                    null,
//                    null,
//                    null,
//                    node.pruneMat,
//                    false,
//                    new EPICSettings(),
//                    null
//            );
//
//            FullAStarNode rootNode = node.stateTree.rootNode();
//
//            int blankConf[] = new int[stateNumPos];//set up root node UB
//            Arrays.fill(blankConf,-1);
//            rootNode.UBConf = blankConf;
//            updateUB(rootNode, node );
//            node.stateUB = rootNode.UB;
//            node.stateTree.initQueue(rootNode);//allocate queue and add root node
//        }
//        else {//no confs available for this state!
//            node.stateTree = null;
//            node.stateUB = Double.POSITIVE_INFINITY;
//        }
//    }
//
//
//
//    @Override
//    public EWAKStarNode rootNode() {
//        int[] conf = new int[numTreeLevels];
//        Arrays.fill(conf,-1);//indicates sequence not assigned
//
//        EWAKStarNode root = new EWAKStarNode(conf, pmat);
//        root.setScore( boundLME(root) );
//        return root;
//    }
//
//
//    @Override
//    public boolean canPruneNode(EWAKStarNode seqNode){
//        //we have a cap on the number of mutations...prune if exceeded or not matched.
//
//        seqNode.getNodeAssignments();
//        if(numMutable<numTreeLevels) {
//            int mutCount = 0;
//            int assignments[] = seqNode.getNodeAssignments();
//
//            int seqListIndex = 0;
//            for (int level = 0; level < numTreeLevels; level++) {
//                if (assignments[level] >= 0) {//AA type at level is assigned
//                    if(confSpace.positions.get(level).resTypes.size()>1) { //assume that position is wild-type if there is only one residue allowed
//                        if (!AATypeOptions.get(level).get(assignments[level]).equalsIgnoreCase(wtSeqList[seqListIndex])){//and is different from wtSeq
//                            mutCount++;
//                            seqListIndex++;
//                        }
//                    }
//                } else {
//                    break;
//                }
//            }
//
//            if (mutCount > numMutable)
//                return true;
//            else if (useExact && seqNode.isFullyDefined() && mutCount != 0 && mutCount != numMutable)
//                return true;
//        }
//        return false;
//
//    }
//
//
//    @Override
//    public boolean isFullyAssigned(EWAKStarNode seqNode) {
//        //This checks if the node is returnable
//        //So it must be fully processed (state GMECs found, not just fully defined sequence)
//
//        if( ! seqNode.isFullyDefined() )//sequence not fully defined
//            return false;
//
//        if(seqNode.stateTree!=null){
//            FullAStarNode bestNodeForState = seqNode.stateTree.getQueue().peek();
//
//            if( ! bestNodeForState.isFullyDefined() )//State GMEC calculation not done
//                return false;
//        }
//
//
//        //if we get here, all state GMECs are calculated.  Node is fully processed
//        return true;
//    }
//
//
//    @Override
//    public ConfSearch.ScoredConf outputNode(EWAKStarNode myNode){
//        //Let's print more info when outputting a node
//        String childSeq = seqAsString(myNode.getNodeAssignments());
//        if (childSeq.equals(NewEWAKStarDoer.toWildTypeEWAKStarString(wtSeq))){
//            System.out.println("Setting WT score as minimized energy...");
//            int conf[] = myNode.stateTree.getQueue().peek().getNodeAssignments();
//            RCTuple rct = new RCTuple(conf);
//            wtMinimizedEnergy = confECalc.calcEnergy(rct).energy;
//        }
//        //printBestSeqInfo(myNode);
//        return new ConfSearch.ScoredConf(myNode.getNodeAssignments(), myNode.getScore());
//    }
//
//
//    public String seqAsString(int[] seqNodeAssignments){
//        //Given the integer representation of sequence used in EWAKStarNode,
//        //create a string representation of the sequence
//        String seq = "";
//        for(int level=0; level<numTreeLevels; level++)
//            if(AATypeOptions.get(level).size()>1) {
//                seq = seq + AATypeOptions.get(level).get(seqNodeAssignments[level]) + "_";
//            }
//
//        return seq;
//    }
//
//
//    private double boundLME(EWAKStarNode seqNode){
//        //Lower-bound func over the sequence space defined by this node
//        if(seqNode.isFullyDefined())//fully-defined sequence
//            return calcLBConfTrees(seqNode);
//        else
//            return calcLBPartialSeq(seqNode);
//    }
//
//
//    private int[] AAOptions(int pos, int assignment){
//        //which of the AA option indices for this mutable position are allowed
//        //given this assignment in nodeAssignments?
//
//        if(assignment==-1){//all options allowed
//            int numOptions = AATypeOptions.get(pos).size();
//            int[] ans = new int[numOptions];
//            for(int option=0; option<numOptions; option++)
//                ans[option] = option;
//
//            return ans;
//        }
//        else//just the assigned option
//            return new int[] {assignment};
//    }
//
//
//    private double calcLBPartialSeq(EWAKStarNode seqNode){
//
//        int partialSeq[] = seqNode.getNodeAssignments();
//
//        double lowerBound = 0;
//        //first terms for mutable residues
//        for(int i=0; i<numTreeLevels; i++){
//            double resE = Double.POSITIVE_INFINITY;
//
//            for( int curAA : AAOptions(i, partialSeq[i]) ){
//                double AAE = 0;
//
//                //residue number i converted to this state's flexible position numbering
//                int statePosNum = mutablePosNums.get(i);
//
//                double stateAAE = Double.POSITIVE_INFINITY;
//
//                PruningMatrix pruneMat = seqNode.pruneMat;
//
//                ArrayList<Integer> rotList = unprunedRCsAtAA(pruneMat, i, statePosNum, curAA);
//
//                for(int rot : rotList){
//
//                    double rotE = emat.getOneBody(statePosNum, rot);
//
//                    for(int pos2=0; pos2<stateNumPos; pos2++){//all non-mut; seq only if < this one
//                        if( (!mutablePosNums.contains(pos2)) || pos2<statePosNum){
//
//                            double bestInteraction = Double.POSITIVE_INFINITY;
//                            ArrayList<Integer> rotList2 = pruneMat.unprunedRCsAtPos(pos2);
//
//                            for(int rot2 : rotList2){
//                                //rot2 known to be unpruned
//                                if(!pruneMat.getPairwise(statePosNum, rot, pos2, rot2)){
//
//                                    double pairwiseE = emat.getPairwise(statePosNum, rot, pos2, rot2);
//                                    pairwiseE += higherOrderContrib(pruneMat, statePosNum, rot, pos2, rot2);
//
//                                    bestInteraction = Math.min(bestInteraction,pairwiseE);
//
//                                }
//                            }
//
//                            rotE += bestInteraction;
//
//                        }
//                    }
//
//                    stateAAE = Math.min(stateAAE,rotE);
//
//                }
//
//                if(Double.isInfinite(stateAAE))
//                    //this will occur if the state is impossible (all confs pruned)
//                    AAE = Double.POSITIVE_INFINITY;
//
//                else
//                    AAE += stateAAE;
//
//
//                resE = Math.min(resE,AAE);
//            }
//
//            lowerBound += resE;
//        }
//
//        //now we bound the energy for the other residues for each of the states
//        //(internal energy for that set of residues)
//        lowerBound += emat.getConstTerm();
//
//        double nonMutBound = boundStateNonMutE(seqNode);
//        //handle this like AAE above
//        if(Double.isInfinite(nonMutBound)){
//            //this will occur if the state is impossible (all confs pruned)
//            return Double.POSITIVE_INFINITY;
//        }
//        else
//            lowerBound += nonMutBound;
//
//        return lowerBound;
//    }
//
//
//    private double higherOrderContrib(PruningMatrix pruneMat,
//                              int pos1, int rc1, int pos2, int rc2){
//        //higher-order contribution for a given RC pair in a given state,
//        //when scoring a partial conf
//
//        HigherTupleFinder<Double> htf = emat.getHigherOrderTerms(pos1,rc1,pos2,rc2);
//
//        if(htf==null)
//            return 0;//no higher-order interactions
//        else{
//            RCTuple curPair = new RCTuple(pos1,rc1,pos2,rc2);
//            return higherOrderContrib(pruneMat, htf, curPair);
//        }
//    }
//
//    private boolean posComesBefore(int pos1, int pos2){
//        //Does pos1 come "before" pos2?  (Ordering: non-mut pos, then mut pos, in ascending order)
//        //These are flexible positions for the given state
//        boolean isMut1 = mutablePosNums.contains(pos1);
//        boolean isMut2 = mutablePosNums.contains(pos2);
//
//        if(isMut1 && !isMut2)
//            return false;
//        if(isMut2 && !isMut1)
//            return true;
//
//        //ok they are in the same category (mut or not)
//        return pos1<pos2;
//    }
//
//
//    private double higherOrderContrib(PruningMatrix pruneMat, HigherTupleFinder<Double> htf,
//                              RCTuple startingTuple){
//        //recursive function to get bound on higher-than-pairwise terms
//        //this is the contribution to the bound due to higher-order interactions
//        //of the RC tuple startingTuple (corresponding to htf)
//
//        double contrib = 0;
//
//        //to avoid double-counting, we are just counting interactions of starting tuple
//        //with residues before the "earliest" one (startingLevel) in startingTuple
//        //"earliest" means lowest-numbered, except non-mutating res come before mutating
//        int startingLevel = startingTuple.pos.get( startingTuple.pos.size()-1 );
//
//        for(int iPos : htf.getInteractingPos() ){//position has higher-order interaction with tup
//            if(posComesBefore(iPos,startingLevel)) {//interaction in right order
//                //(want to avoid double-counting)
//
//                double levelBestE = Double.POSITIVE_INFINITY;//best value of contribution
//                //from tup-iPos interaction
//
//                ArrayList<Integer> allowedRCs = pruneMat.unprunedRCsAtPos(iPos);
//
//                for( int rc : allowedRCs ){
//
//                    RCTuple augTuple = startingTuple.addRC(iPos, rc);
//
//                    if( ! pruneMat.isPruned(augTuple) ){
//
//                        double interactionE = htf.getInteraction(iPos, rc);
//
//                        //see if need to go up to highers order again...
//                        HigherTupleFinder<Double> htf2 = htf.getHigherInteractions(iPos, rc);
//                        if(htf2!=null){
//                            interactionE += higherOrderContrib(pruneMat,htf2,augTuple);
//                        }
//
//                        //besides that only residues in definedTuple or levels below level2
//                        levelBestE = Math.min(levelBestE,interactionE);
//
//                    }
//                }
//
//                contrib += levelBestE;//add up contributions from different interacting positions iPos
//            }
//        }
//
//        return contrib;
//    }
//
//
//
//
//    private ArrayList<Integer> unprunedRCsAtAA(PruningMatrix pruneMat,
//                                               int mutablePos, int statePosNum, int curAA){
//        //List the RCs of the given position in the given state
//        //that come with the indicated AA type (indexed in AATypeOptions)
//
//        ArrayList<Integer> unprunedRCs = pruneMat.unprunedRCsAtPos(statePosNum);
//        ArrayList<Integer> ans = new ArrayList<>();
//
//        String curAAType = AATypeOptions.get(mutablePos).get(curAA);
//
//        for(int rc : unprunedRCs){
//            String rcAAType = confSpace.positions.get(statePosNum).resConfs.get(rc).template.name;
//
//            if(rcAAType.equalsIgnoreCase(curAAType))//right type
//                ans.add(rc);
//        }
//
//        return ans;
//    }
//
//
//
//
//    private double calcLBConfTrees(EWAKStarNode seqNode){
//        //here the sequence is fully defined
//        //so we can bound func solely based on lower and upper bounds (depending on func coefficient sign)
//        //of the GMECs for each state, which can be derived from the front node of each state's ConfTree
//        double ans = 0;
//
//
//        if(seqNode.stateTree==null)//state and sequence impossible
//            return Double.POSITIVE_INFINITY;
//
//        ans += seqNode.stateTree.getQueue().peek().getScore();
//
//
//        //shell-shell energies may differ between states!
//        //THIS DOUBLE-COUNTS BC SCORES ALREADY INCLUDE CONST TERM
//        //ans += func.coeffs[state]*getEnergyMatrix(state).getConstTerm();
//
//
//        return ans;
//    }
//
//
//    private double boundStateNonMutE(EWAKStarNode seqNode){
//        //get a quick lower or upper bound (as indicated) for the energy of the given state's
//        //non-mutable residues (their intra+shell energies + pairwise between them)
//        //use pruning information from seqNode
//        double ans = 0;
//
//        PruningMatrix pruneMat = seqNode.pruneMat;
//
//        for(int pos=0; pos<stateNumPos; pos++){
//            if((!mutablePosNums.contains(pos))){
//
//                double resE = Double.POSITIVE_INFINITY;
//
//                ArrayList<Integer> rotList = pruneMat.unprunedRCsAtPos(pos);
//
//
//                for(int rot : rotList){
//                    //make sure rot isn't pruned
//                    if(!pruneMat.getOneBody(pos, rot)){
//
//                        double rotE = emat.getOneBody(pos, rot);
//
//                        for(int pos2=0; pos2<stateNumPos; pos2++){//all non-mut; seq only if < this one
//                            if( (!mutablePosNums.contains(pos2)) && pos2<pos){
//
//                                double bestInteraction = Double.POSITIVE_INFINITY;
//                                ArrayList<Integer> rotList2 = pruneMat.unprunedRCsAtPos(pos2);
//
//                                for(int rot2 : rotList2){
//                                    if(!pruneMat.getPairwise(pos, rot, pos2, rot2)){
//                                        double pairwiseE = emat.getPairwise(pos,rot,pos2,rot2);
//
//                                        pairwiseE += higherOrderContrib(pruneMat, pos, rot, pos2, rot2);
//
//                                        bestInteraction = Math.min(bestInteraction,pairwiseE);
//
//                                    }
//                                }
//
//                                rotE += bestInteraction;
//                            }
//                        }
//
//                        resE = Math.min(resE,rotE);
//
//                    }
//                }
//
//                ans += resE;
//            }
//        }
//
//        return ans;
//    }
//
//
//    private boolean updateUB(FullAStarNode expNode, EWAKStarNode seqNode){
//        //Get an upper-bound on the node by a little FASTER run, generating UBConf
//        //store UBConf and UB in expNode
//        //expNode is in seqNode.stateTrees[state]
//        //we'll start with the starting conf (likely from a parent) if provided
//        //return true unless no valid conf is possible...then false
//
//
//        int assignments[] = expNode.getNodeAssignments();
//        ArrayList<ArrayList<Integer>> allowedRCs = new ArrayList<>();
//
//        for(int pos=0; pos<stateNumPos; pos++){
//
//            ArrayList<Integer> posOptions = new ArrayList<>();
//
//            if(assignments[pos]==-1){
//                posOptions = seqNode.pruneMat.unprunedRCsAtPos(pos);
//                if(posOptions.isEmpty())//no options at this position
//                    return false;
//            }
//            else//assigned
//                posOptions.add(assignments[pos]);
//
//            allowedRCs.add(posOptions);
//        }
//
//
//        int startingConf[] = expNode.UBConf;
//        int[] UBConf = startingConf.clone();
//        //ok first get rid of anything in startingConf not in expNode's conf space,
//        //replace with the lowest-intra+shell-E conf
//        for(int level=0; level<stateNumPos; level++){
//
//            if( ! allowedRCs.get(level).contains(startingConf[level]) ){
//                //if( ! levelOptions.get(level).get(expNode.conf[level]).contains( startingConf[level] ) ){
//
//                double bestISE = Double.POSITIVE_INFINITY;
//                int bestRC = allowedRCs.get(level).get(0);
//                for(int rc : allowedRCs.get(level) ){
//                    double ise = emat.getOneBody(level, rc);
//                    if( ise < bestISE){
//                        bestISE = ise;
//                        bestRC = rc;
//                    }
//                }
//
//                UBConf[level] = bestRC;
//            }
//        }
//
//
//        double curE = emat.confE(UBConf);
//        boolean done = false;
//
//        while(!done){
//
//            done = true;
//
//            for(int level=0; level<stateNumPos; level++){
//
//                int testConf[] = UBConf.clone();
//
//                for(int rc : allowedRCs.get(level) ){
//                    testConf[level] = rc;
//
//                    //if(!canPrune(testConf)){//pruned conf unlikely to be good UB
//                    //would only prune if using pruned pair flags in A*
//
//                    double testE = emat.confE(testConf);
//                    if( testE < curE){
//                        curE = testE;
//                        UBConf[level] = rc;
//                        done = false;
//                    }
//                    //}
//                }
//            }
//        }
//
//
//        expNode.UBConf = UBConf;
//        expNode.UB = emat.confE(UBConf);
//
//        return true;
//    }
//
//}