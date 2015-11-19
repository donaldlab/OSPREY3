/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar;

import edu.duke.cs.osprey.confspace.ConfSpaceSuper;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.SearchProblemSuper;
import edu.duke.cs.osprey.confspace.SuperRCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author hmn5
 */
//ConfTree Class but with super-RCs
public class ConfTreeSuper extends AStarTree {
    //This implementation of an A* tree is intended for conformational search
    //AStarNode.nextAssignment is an array of length numPos; each position
    //stores the assigned RC, or -1 to indicate an unassigned position
    //this class supports both "traditional" (static, simple heuristic) A*
    //and improvements like dynamic A*
    //we may also want to allow other negative indices, to indicate partially assigned RCs

    int numPos;
    EnergyMatrix emat;

    ArrayList<ArrayList<Integer>> unprunedRCsAtPos = new ArrayList<>();
    //get from searchSpace when initializing!
    //These are lists of residue-specific RC numbers for the unpruned RCs at each residue

    //ADVANCED SCORING METHODS: TO CHANGE LATER (EPIC, MPLP, etc.)
    boolean traditionalScore = true;
    boolean mplpScore = true;
    boolean useRefinement = true;//refine nodes (might want EPIC, MPLP, or something else)
    boolean useEpic;

    boolean useDynamicAStar = true;

    EPICMatrix epicMat = null;//to use in refinement
    ConfSpaceSuper confSpace = null;//conf space to use with epicMat if we're doing EPIC minimization w/ SAPE
    boolean minPartialConfs = false;//whether to minimize partially defined confs with EPIC, or just fully defined

    // MPLP object for node refinement.
    public Mplp mplpMinimizer;

    public ConfTreeSuper(SearchProblemSuper sp) {
        init(sp, sp.pruneMat, sp.useEPIC);
    }

    public ConfTreeSuper(SearchProblemSuper sp, PruningMatrix pruneMat, boolean useEPIC) {
        //Conf search over RC's in sp that are unpruned in pruneMat
        init(sp, pruneMat, useEPIC);
    }

    private void init(SearchProblemSuper sp, PruningMatrix pruneMat, boolean useEPIC) {
        numPos = sp.confSpaceSuper.numPos;

        //see which RCs are unpruned and thus available for consideration
        for (int pos = 0; pos < numPos; pos++) {
            unprunedRCsAtPos.add(pruneMat.unprunedRCsAtPos(pos));
        }
        this.useEpic = useEPIC;

        //get the appropriate energy matrix to use in this A* search
        if (sp.useTupExpForSearch) {
            emat = sp.tupExpEMat;
        } else {
            emat = sp.emat;

            if (useEPIC) {//include EPIC in the search
                useRefinement = true;
                epicMat = sp.epicMat;
                confSpace = sp.confSpaceSuper;
                minPartialConfs = sp.epicSettings.minPartialConfs;
            }
        }
        // Initialize MPLP
        if (mplpScore) {
            useRefinement = true;
            mplpMinimizer = new Mplp(numPos, unprunedRCsAtPos, emat);
        }
    }

    @Override
    public ArrayList<AStarNode> getChildren(AStarNode curNode) {

        if (isFullyAssigned(curNode)) {
            throw new RuntimeException("ERROR: Can't expand a fully assigned A* node");
        }

        ArrayList<AStarNode> ans = new ArrayList<>();
        int nextLevel = nextLevelToExpand(curNode.nodeAssignments);

        for (int rc : unprunedRCsAtPos.get(nextLevel)) {
            int[] childConf = curNode.nodeAssignments.clone();
            childConf[nextLevel] = rc;
            AStarNode childNode = new AStarNode(childConf, scoreConf(childConf), useRefinement);
            ans.add(childNode);
        }

        return ans;
    }

    @Override
    public AStarNode rootNode() {
        //no residues assigned, so all -1's
        int[] conf = new int[numPos];
        Arrays.fill(conf, -1);

        AStarNode root = new AStarNode(conf, scoreConf(conf), useRefinement);
        return root;
    }

    @Override
    public boolean isFullyAssigned(AStarNode node) {
        for (int rc : node.nodeAssignments) {
            if (rc < 0)//not fully assigned
            {
                return false;
            }
        }

        return true;
    }

    //operations supporting special features like dynamic A*
    /**
     * given a partially defined conformation, what level should be expanded
     * next?
     *
     * @param partialConf Partially defined rotamer conformation
     * @return
     */
    public int nextLevelToExpand(int[] partialConf) {

        if (useDynamicAStar) {

            int bestLevel = -1;
            double bestLevelScore = Double.NEGATIVE_INFINITY;

            for (int level = 0; level < numPos; level++) {
                if (partialConf[level] < 0) {//position isn't already all expanded

                    double levelScore = scoreExpansionLevel(level, partialConf);

                    if (levelScore > bestLevelScore) {//higher score is better
                        bestLevelScore = levelScore;
                        bestLevel = level;
                    }
                }
            }

            if (bestLevel == -1) {
                throw new RuntimeException("ERROR: No next expansion level found for dynamic A*");
            }

            return bestLevel;
        } else {//static ordering.  
            //Let's only support the traditional ordering since dynamic will beat static for improved orderings.
            for (int level = 0; level < numPos; level++) {
                if (partialConf[level] < 0) {
                    return level;
                }
            }

            throw new RuntimeException("ERROR: Can't find next expansion level for fully defined conformation");
        }

    }

    /**
     * Score expansion at the indicated level for the given partial conformation
     * for use in dynamic A*. Higher score is better.
     *
     * @param level
     * @param partialConf
     * @return
     */
    double scoreExpansionLevel(int level, int[] partialConf) {

        //best performing score is just 1/(sum of reciprocals of score rises for child nodes)
        double parentScore = scoreConf(partialConf);
        int[] expandedConf = partialConf.clone();

        double reciprocalSum = 0;

        for (int rc : unprunedRCsAtPos.get(level)) {
            expandedConf[level] = rc;
            double childScore = scoreConf(expandedConf);

            reciprocalSum += 1.0 / (childScore - parentScore);
        }

        double score = 1. / reciprocalSum;

        return score;
    }

    double scoreConf(int[] partialConf) {
        if (traditionalScore) {
            SuperRCTuple definedTuple = new SuperRCTuple(partialConf);

            double score = emat.getConstTerm() + emat.getInternalEnergy(definedTuple);//"g-score"

            //score works by breaking up the full energy into the energy of the defined set of residues ("g-score"),
            //plus contributions associated with each of the undefined res ("h-score")
            for (int level = 0; level < numPos; level++) {
                if (partialConf[level] < 0) {//level not fully defined

                    double resContribLB = Double.POSITIVE_INFINITY;//lower bound on contribution of this residue
                    //resContribLB will be the minimum_{rc} of the lower bound assuming rc assigned to this level

                    for (int rc : unprunedRCsAtPos.get(level)) {
                        resContribLB = Math.min(resContribLB, RCContributionLB(level, rc, definedTuple, partialConf));
                    }

                    score += resContribLB;
                }
            }

            return score;
        } else {

            //other possibilities include MPLP, etc.
            //But I think these are better used as refinements
            //we may even want multiple-level refinement
            throw new RuntimeException("Advanced A* scoring methods not implemented yet!");
        }
    }

    /**
     * Provide a lower bound on what the given rc at the given level can
     * contribute to the energy assume partialConf and definedTuple
     *
     * @param level
     * @param rc
     * @param definedTuple
     * @param partialConf
     * @return
     */
    double RCContributionLB(int level, int rc, SuperRCTuple definedTuple, int[] partialConf) {
        double rcContrib = emat.getOneBody(level, rc);

        //for this kind of lower bound, we need to split up the energy into the defined-tuple energy
        //plus "contributions" for each undefined residue
        //so we'll say the "contribution" consists of any interactions that include that residue
        //but do not include higher-numbered undefined residues
        for (int level2 = 0; level2 < numPos; level2++) {

            if (partialConf[level2] >= 0 || level2 < level) {//lower-numbered or defined residues

                double levelBestE = Double.POSITIVE_INFINITY;//best pairwise energy
                ArrayList<Integer> allowedRCs = allowedRCsAtLevel(level2, partialConf);

                for (int rc2 : allowedRCs) {

                    double interactionE = emat.getPairwise(level, rc, level2, rc2);

                    double higherLB = higherOrderContribLB(partialConf, level, rc, level2, rc2);
                    //add higher-order terms that involve rc, rc2, and parts of partialConf

                    interactionE += higherLB;

                    //besides that only residues in definedTuple or levels below level2
                    levelBestE = Math.min(levelBestE, interactionE);
                }

                rcContrib += levelBestE;
            }
        }

        return rcContrib;
    }

    ArrayList<Integer> allowedRCsAtLevel(int level, int[] partialConf) {
        //What RCs are allowed at the specified level (i.e., position num) in the given partial conf?
        ArrayList<Integer> allowedRCs;

        if (partialConf[level] == -1)//position undefined: consider all RCs
        {
            allowedRCs = unprunedRCsAtPos.get(level);
        } else if (partialConf[level] >= 0) {
            allowedRCs = new ArrayList<>();
            allowedRCs.add(partialConf[level]);
        } else {
            throw new UnsupportedOperationException("ERROR: Partially assigned position not yet supported in A*");
        }

        return allowedRCs;
    }

    double higherOrderContribLB(int[] partialConf, int pos1, int rc1, int pos2, int rc2) {
        //higher-order contribution for a given RC pair, when scoring a partial conf

        HigherTupleFinder<Double> htf = emat.getHigherOrderTerms(pos1, rc1, pos2, rc2);

        if (htf == null) {
            return 0;//no higher-order interactions
        } else {
            return higherOrderContribLB(partialConf, htf, pos2);
        }
    }

    /**
     * recursive function to get lower bound on higher-than-pairwise terms this
     * is the contribution to the lower bound due to higher-order interactions
     * of the RC tuple corresponding to htf with "lower-numbered" residues
     * (numbering as in scoreConf: these are residues that are fully defined in
     * partialConf, or are actually numbered <level2) @
     *
     *
     * param partialConf @param htf @param level2
     * @r
     *
     * e
     * t
     * urn
     */
    double higherOrderContribLB(int[] partialConf, HigherTupleFinder<Double> htf, int level2) {

        double contrib = 0;

        for (int iPos : htf.getInteractingPos()) {//position has higher-order interaction with tup
            if (posComesBefore(iPos, level2, partialConf)) {//interaction in right order
                //(want to avoid double-counting)

                double levelBestE = Double.POSITIVE_INFINITY;//best value of contribution
                //from tup-iPos interaction
                ArrayList<Integer> allowedRCs = allowedRCsAtLevel(iPos, partialConf);

                for (int rc : allowedRCs) {

                    double interactionE = htf.getInteraction(iPos, rc);

                    //see if need to go up to highers order again...
                    HigherTupleFinder htf2 = htf.getHigherInteractions(iPos, rc);
                    if (htf2 != null) {
                        interactionE += higherOrderContribLB(partialConf, htf2, iPos);
                    }

                    //besides that only residues in definedTuple or levels below level2
                    levelBestE = Math.min(levelBestE, interactionE);
                }

                contrib += levelBestE;//add up contributions from different interacting positions iPos
            }
        }

        return contrib;
    }

    /**
     * posComesBefore: for purposes of contributions to traditional conf score,
     * we go through defined and then through undefined positions (in
     * partialConf); within each of these groups we go in order of position
     * number
     *
     * @param pos1
     * @param pos2
     * @param partialConf
     * @return
     */
    private boolean posComesBefore(int pos1, int pos2, int partialConf[]) {

        if (partialConf[pos2] >= 0) {//pos2 defined
            return (pos1 < pos2 && partialConf[pos1] >= 0);//pos1 must be defined to come before pos2
        } else//pos1 comes before pos2 if it's defined, or if pos1<pos2
        {
            return (pos1 < pos2 || partialConf[pos1] >= 0);
        }
    }

    @Override
    void refineScore(AStarNode node) {

        if (useEpic && epicMat == null) {
            throw new UnsupportedOperationException("ERROR: Trying to call refinement w/o EPIC matrix");
        }
        // Refine node with MPLP.
        if (this.mplpScore) {
            node.score = this.mplpMinimizer.optimizeEMPLP(node.nodeAssignments, 100);
        }

        if (useEpic && (minPartialConfs || isFullyAssigned(node))) {
            node.score += epicMat.minContE(node.nodeAssignments);
        }

        node.scoreNeedsRefinement = false;
    }

    public void setMPLPForInteractionEnergy(List<EnergyMatrix> boundResNumToUnboundEmat, List<Integer> boundResNumToUnboundResNum,
            List<Boolean> boundresNumToIsMutableStrand, boolean[][] belongToSameStrand) {
        this.mplpMinimizer.setCrossTermInteractionGraph(boundResNumToUnboundEmat, boundResNumToUnboundResNum, boundresNumToIsMutableStrand, belongToSameStrand);
        this.mplpMinimizer.createOnlyPairwiseMat(this.emat, boundResNumToUnboundEmat, boundResNumToUnboundResNum, boundresNumToIsMutableStrand, belongToSameStrand);
    }

    //this function computes the minimum over all full conf E's consistent with partialConf
    //for debugging only of course
    double exhaustiveScore(int[] partialConf) {
        for (int pos = 0; pos < partialConf.length; pos++) {
            if (partialConf[pos] == -1) {
                //recurse to get all options
                double score = Double.POSITIVE_INFINITY;
                for (int rc : allowedRCsAtLevel(pos, partialConf)) {
                    int partialConf2[] = partialConf.clone();
                    partialConf2[pos] = rc;
                    score = Math.min(score, exhaustiveScore(partialConf2));
                }
                return score;
            }
        }
        //if we get here, conf fully defined
        return emat.getInternalEnergy(new SuperRCTuple(partialConf));
    }

    public double energyNextConf() {
        return this.scoreConf(this.nextConf());
    }

    /**
     * Compute a polynomial time, upper bound on the energy of a partial conformation in the conformation tree using MPLP.  
     * @param conf: 
     * @return
     */
    public double mplpUpperBound(int[] conf) {
        return this.mplpMinimizer.optimizeEMPLP(conf, 100);
    }
}
