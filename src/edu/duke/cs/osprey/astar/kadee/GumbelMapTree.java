/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.kadee;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.AStarTree;
import edu.duke.cs.osprey.astar.Mplp;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.partitionfunctionbounds.GumbelDistribution;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author hmn5
 */
public class GumbelMapTree extends AStarTree {

    int numPos;
    EnergyMatrix emat;
    //HMN: Added Pruning Mat for Pairs Pruning
    PruningMatrix pruneMat;

    ArrayList<ArrayList<Integer>> unprunedRCsAtPos = new ArrayList<>();
    //get from searchSpace when initializing!
    //These are lists of residue-specific RC numbers for the unpruned RCs at each residue

    //ADVANCED SCORING METHODS: TO CHANGE LATER (EPIC, MPLP, etc.)
    public boolean traditionalScore = true;
    boolean useRefinement = true;
    public boolean mplpScore = false;

    //MPLP object for node refinement
    public Mplp mplpMinimizer;

    boolean useDynamicAStar = true;

    ConfSpace confSpace = null;//conf space to use with epicMat if we're doing EPIC minimization w/ SAPE
    boolean minPartialConfs = false;//whether to minimize partially defined confs with EPIC, or just fully defined

    public double epsilon = 1.0;

    int[] currentBestFeasibleSolution = null;
    public double currentBestFeasibleScore = Double.POSITIVE_INFINITY;
    public double currentBestLPRelaxScore = Double.POSITIVE_INFINITY;

    public double lowerBoundLogZ = Double.NEGATIVE_INFINITY;
    public double upperBoundLogZ = Double.POSITIVE_INFINITY;
    public int numNodesEpsilon = -1;

    private Random randomGenerator;
    final double constRT = PoissonBoltzmannEnergy.constRT;

    boolean verbose = false;

    
    public GumbelMapTree(SearchProblem sp) {
        init(sp, sp.pruneMat);
    }

    public GumbelMapTree(EnergyMatrix aEmat, PruningMatrix aPruneMat) {
        this.numPos = aEmat.numPos();
        this.emat = aEmat;
        this.pruneMat = aPruneMat;

        //see which RCs are unpruned and thus available for consideration
        for (int pos = 0; pos < numPos; pos++) {
            unprunedRCsAtPos.add(aPruneMat.unprunedRCsAtPos(pos));
        }

        if (mplpScore) {
            mplpMinimizer = new Mplp(numPos, emat, pruneMat);
        }
        randomGenerator = new Random();
    }
    
    
    private void init(SearchProblem sp, PruningMatrix aPruneMat) {
        numPos = sp.confSpace.numPos;

        //see which RCs are unpruned and thus available for consideration
        for (int pos = 0; pos < numPos; pos++) {
            unprunedRCsAtPos.add(aPruneMat.unprunedRCsAtPos(pos));
        }
        this.pruneMat = aPruneMat;
        //get the appropriate energy matrix to use in this A* search
        if (sp.useTupExpForSearch) {
            emat = sp.tupExpEMat;
        } else {
            emat = sp.emat;
            if (sp.useEPIC) {//include EPIC in the search
                throw new RuntimeException("ERROR: GumbelMapTree Does Not Support EPIC");
            }
        }
        //Initialize MPLP (not compatible with EPIC)
        if (mplpScore) {
            //mplpMinimizer = new Mplp(numPos, emat, pruneMat);
            mplpMinimizer = new Mplp(this.numPos, emat, pruneMat);
            if (sp.useTupExpForSearch) {
                throw new RuntimeException("ERROR: MPLP Does Not Yet work with LUTE");
            }
        }
        randomGenerator = new Random();
    }

    @Override
    public ArrayList<AStarNode> getChildren(AStarNode curNode) {

        if (isFullyAssigned(curNode)) {
            throw new RuntimeException("ERROR: Can't expand a fully assigned A* node");
        }

        if (curNode.getScore() == Double.POSITIVE_INFINITY)//node impossible, so no children
        {
            return new ArrayList<>();
        }

        ArrayList<AStarNode> ans = new ArrayList<>();
        int nextLevel = -1;
        if (!nodeAssigned(curNode)) {
            nextLevel = nextLevelToExpand(curNode.getNodeAssignments(), curNode.getScore(), curNode.perturbation);
        } else {
            return ans;
        }
        for (int rc : unprunedRCsAtPos.get(nextLevel)) {
            int[] childConf = curNode.getNodeAssignments().clone();
            childConf[nextLevel] = rc;

            double logSearchProblemSize = computeLogSearchSize(childConf);

            //Pass perturbation down
            if (curNode.getNodeAssignments()[nextLevel] == rc) {
                double score = scoreConfWithPert(childConf, curNode.perturbation);
                AStarNode childNode = new AStarNode(childConf, score, curNode.feasibleSolution, curNode.perturbation, useRefinement);
                if (!canPruneNodeGumbel(childNode)){
                    ans.add(childNode);
                }
                updateCurrentBestFeasibleSolution(childNode);
            } else {
                double gumbelPert = GumbelDistribution.sampleTruncated(-GumbelDistribution.gamma + logSearchProblemSize, curNode.perturbation);
                double score = scoreConfWithPert(childConf, gumbelPert);
                AStarNode childNode = new AStarNode(childConf, score, getFeasibleSolution(childConf), gumbelPert, useRefinement);
                if (!canPruneNodeGumbel(childNode)){
                    ans.add(childNode);
                }
                updateCurrentBestFeasibleSolution(childNode);
            }
        }
        return ans;
    }

    double scoreConf(int[] partialConf) {
        if (confAssigned(partialConf)) {
            return emat.getConstTerm() + emat.getInternalEnergy(new RCTuple(partialConf));
        }

        if (traditionalScore) {
            RCTuple definedTuple = new RCTuple(partialConf);

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
        } else if (mplpScore) {
            double score = emat.getConstTerm() + this.mplpMinimizer.optimizeMPLP(partialConf, 100);
//            double score = emat.getConstTerm() + this.mplpMinimizer.optimizeEMPLP(partialConf, numPos);
//            double scoreTraditional = scoreConfTraditional(partialConf);

            return score;

        } else {
            throw new RuntimeException("Advanced A* scoring methods not implemented yet!");
        }

    }

    public double scoreConfTraditional(int[] partialConf) {
        RCTuple definedTuple = new RCTuple(partialConf);

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
    }

    public double scoreConfTraditionalWithPert(int[] partialConf, double pert) {
        RCTuple definedTuple = new RCTuple(partialConf);

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

        return score - this.constRT * pert;
    }

    //operations supporting special features like dynamic A*
    public int nextLevelToExpand(int[] partialConf, double curLevelScore, double parentNoise) {
        //given a partially defined conformation, what level should be expanded next?

        if (useDynamicAStar) {

            int bestLevel = -1;
            double bestLevelScore = Double.NEGATIVE_INFINITY;

            for (int level = 0; level < numPos; level++) {
                if (partialConf[level] < 0) {//position isn't already all expanded

                    double levelScore = scoreExpansionLevel(level, partialConf, curLevelScore, parentNoise);

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

    double scoreExpansionLevel(int level, int[] partialConf, double curLevelScore, double parentNoise) {
        //Score expansion at the indicated level for the given partial conformation
        //for use in dynamic A*.  Higher score is better.

        //best performing score is just 1/(sum of reciprocals of score rises for child nodes)
        double parentScore = scoreConfTraditionalWithPert(partialConf, parentNoise);
        int[] expandedConf = partialConf.clone();

        double reciprocalSum = 0;

        for (int rc : unprunedRCsAtPos.get(level)) {
            expandedConf[level] = rc;
            double logConfSpace = computeLogSearchSize(expandedConf);
            double gumbelNoise = GumbelDistribution.sampleTruncated(-GumbelDistribution.gamma + logConfSpace, parentNoise);
//            double childScore = scoreConf(expandedConf);
            double childScore = scoreConfTraditionalWithPert(expandedConf, gumbelNoise);
            reciprocalSum += 1.0 / (childScore - parentScore);
        }

        double score = 1. / reciprocalSum;

        return score;
    }

    double RCContributionLB(int level, int rc, RCTuple definedTuple, int[] partialConf) {
        //Provide a lower bound on what the given rc at the given level can contribute to the energy
        //assume partialConf and definedTuple

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

    double higherOrderContribLB(int[] partialConf, HigherTupleFinder<Double> htf, int level2) {
        //recursive function to get lower bound on higher-than-pairwise terms
        //this is the contribution to the lower bound due to higher-order interactions
        //of the RC tuple corresponding to htf with "lower-numbered" residues (numbering as in scoreConf:
        //these are residues that are fully defined in partialConf, or are actually numbered <level2)

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

    private boolean posComesBefore(int pos1, int pos2, int partialConf[]) {
        //for purposes of contributions to traditional conf score, 
        //we go through defined and then through undefined positions (in partialConf);
        //within each of these groups we go in order of position number
        if (partialConf[pos2] >= 0) {//pos2 defined
            return (pos1 < pos2 && partialConf[pos1] >= 0);//pos1 must be defined to come before pos2
        } else//pos1 comes before pos2 if it's defined, or if pos1<pos2
        {
            return (pos1 < pos2 || partialConf[pos1] >= 0);
        }
    }

    @Override
    public void refineScore(AStarNode node) {
    //I will use refineScore to update the bounds 
        updateBounds(node);
    }

    @Override
    public boolean isFullyAssigned(AStarNode node) {
        /*
         for (int rc : node.nodeAssignments) {
         if (rc < 0)//not fully assigned
         {
         return false;
         }
         }
         this.currentBestFeasibleScore = node.score;
         this.currentBestFeasibleSolution = node.nodeAssignments;
         return true;
         */
        return this.pq.isEmpty() && (!node.isRoot);
    }

    @Override
    public AStarNode rootNode() {
        //no residues assigned, so all -1's
        int[] conf = new int[numPos];
        Arrays.fill(conf, -1);

        double logSearchProblemSize = computeLogSearchSize(conf);
        double gumbelNoise = GumbelDistribution.sample(-GumbelDistribution.gamma + logSearchProblemSize, 1);

        int[] feasibleSolution = getFeasibleSolution(conf);
        double lowerBound = scoreConfWithPert(conf, gumbelNoise);

        AStarNode root = new AStarNode(conf, lowerBound, useRefinement);

        root.perturbation = gumbelNoise;
        root.isRoot = true;

        this.currentBestFeasibleSolution = feasibleSolution;
        this.currentBestFeasibleScore = scoreConfWithPert(feasibleSolution, gumbelNoise);

        return root;
    }

    /**
     * Computes the log of the size of the search problem under a particular
     * node,
     *
     * @param partialAssignment the assignments of the node
     * @return
     */
    double computeLogSearchSize(int[] partialAssignment) {
        double size = 0.0;
        for (int pos = 0; pos < this.numPos; pos++) {
            if (partialAssignment[pos] == -1) {
                size += Math.log(this.pruneMat.unprunedRCsAtPos(pos).size());
            }
        }
        return size;
    }

    int[] getFeasibleSolution(int[] partialAssignment) {
        return getRandomConformation(partialAssignment);
        //HMN: proposed change to this.mplpMinimizer.feasibleSolution
    }

    int[] getRandomConformation(int[] partialAssignment) {
        int[] randomConf = new int[partialAssignment.length];
        for (int pos = 0; pos < this.numPos; pos++) {
            if (partialAssignment[pos] == -1) {
                ArrayList<Integer> unprunedRCs = pruneMat.unprunedRCsAtPos(pos);
                int randRotamer = unprunedRCs.get(randomGenerator.nextInt(unprunedRCs.size()));
                randomConf[pos] = randRotamer;
            } else {
                randomConf[pos] = partialAssignment[pos];
            }
        }
        return randomConf;
    }

    double scoreConfWithPert(int[] conf, double perturbation) {
        double confScore = scoreConf(conf);
        return confScore - this.constRT * perturbation;
    }

    public void updateBounds(AStarNode node) {
        this.upperBoundLogZ = -getUpperBoundLogZ(node) / this.constRT;
        this.lowerBoundLogZ = -this.currentBestFeasibleScore / this.constRT;
        if (verbose) {
            System.out.println("Upper Bound logZ: " + this.upperBoundLogZ);
            System.out.println("Lower Bound logZ: " + this.lowerBoundLogZ);
            System.out.println("Gap between Bounds: " + (this.upperBoundLogZ - this.lowerBoundLogZ));
        }
        node.setScoreNeedsRefinement(false);
    }

    public boolean canPruneNodeGumbel(AStarNode node) {
        if (node.getScore() + 0.00001 > this.currentBestFeasibleScore) {
            return true;
        }
        return false;
    }

    public void updateCurrentBestFeasibleSolution(AStarNode node) {
        double score = scoreConfWithPert(node.feasibleSolution, node.perturbation);
        if (score < this.currentBestFeasibleScore) {
            this.currentBestFeasibleScore = score;
            this.currentBestFeasibleSolution = node.feasibleSolution;
        }
    }

    public double getUpperBoundLogZ(AStarNode curNode) {
        double minScore = this.currentBestFeasibleScore;
        for (AStarNode node : this.pq) {
            if (node.getScore() < this.currentBestFeasibleScore) {
                minScore = Math.min(minScore, node.getScore());
            }
        }
        //Also check with current node
        minScore = Math.min(minScore, curNode.getScore());
        return minScore;
    }

    boolean nodeAssigned(AStarNode node) {
        for (int pos : node.getNodeAssignments()) {
            if (pos == -1) {
                return false;
            }
        }
        return true;
    }

    boolean confAssigned(int[] conf) {
        for (int pos : conf) {
            if (pos == -1) {
                return false;
            }
        }
        return true;
    }

    void printNode(AStarNode node) {
        for (int rot : node.getNodeAssignments()) {
            if (rot != -1) {
                System.out.print("|||   ");
            }
        }
        System.out.println();
    }
}
