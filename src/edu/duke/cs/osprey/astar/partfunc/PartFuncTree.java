/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.partfunc;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.AStarTree;
import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.Mplp;
import edu.duke.cs.osprey.astar.comets.UpdatedPruningMatrix;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.partitionfunctionbounds.MapPerturbation;
import edu.duke.cs.osprey.partitionfunctionbounds.ReparamMRF;
import edu.duke.cs.osprey.partitionfunctionbounds.SCMF_Clamp;
import edu.duke.cs.osprey.partitionfunctionbounds.TRBPSeq;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ExpFunction;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author hmn5
 */
public class partFuncTree extends AStarTree {

    double expNormalizer = 0.0;

    BigDecimal lbZ;
    BigDecimal ubZ;

    BigDecimal runningSum;
    double runningSumExpNormalizer = 0.0;

    double epsilon = 0.1;

    int numPos;

    EnergyMatrix emat;
    PruningMatrix pruneMat;

    final double eCut = 0.0; //if we want to set an energy cutoff for MRFs

    boolean useMapPert = false;
    int numSamples = 5; //number of sample averages to take the max from;

    ExpFunction ef = new ExpFunction();
    double constRT = PoissonBoltzmannEnergy.constRT;

    boolean useDynamicOrdering = true;
    boolean branchHalfSpace = true;
    
    Mplp mplp;

    boolean verbose = true;

    public int numConfsEnumerated = 0;

    public partFuncTree(SearchProblem sp) {
        init(sp, sp.pruneMat, sp.useEPIC);
    }

    public partFuncTree(SearchProblem sp, PruningMatrix aPruneMat, boolean useEPIC) {
        //Conf search over RC's in sp that are unpruned in pruneMat
        init(sp, aPruneMat, useEPIC);
    }

    //For rigid traditional
    public partFuncTree(EnergyMatrix aEmat, PruningMatrix aPruneMat) {
        this.numPos = aEmat.numPos();
        this.emat = aEmat;
        this.pruneMat = aPruneMat;

        this.lbZ = new BigDecimal("0.0");
        this.ubZ = new BigDecimal("0.0");
        this.runningSum = new BigDecimal("0.0");
        mplp = new Mplp(numPos, aEmat, aPruneMat);
    }

    private void init(SearchProblem sp, PruningMatrix aPruneMat, boolean useEPIC) {
        numPos = sp.confSpace.numPos;
        this.lbZ = new BigDecimal("0.0");
        this.ubZ = new BigDecimal("0.0");
        this.runningSum = new BigDecimal("0.0");

        this.pruneMat = aPruneMat;
        //get the appropriate energy matrix to use in this A* search
        if (sp.useTupExpForSearch) {
            emat = sp.tupExpEMat;
            throw new RuntimeException("partFuncTree does not yet support tupExpEmat");
        } else {
            emat = sp.emat;

            if (useEPIC) {//include EPIC in the search
                throw new RuntimeException("partFuncTree does not yet support EPIC");
            }
        }
        mplp = new Mplp(numPos, emat, this.pruneMat);
    }

    @Override
    public ArrayList<AStarNode> getChildren(AStarNode curNode) {
        partFuncNode node = (partFuncNode) curNode;
        ArrayList<AStarNode> children = new ArrayList<>();

//        subtractFromBounds(node);
        subtractLowerBound(node);
        int[] curAssignments = node.getNodeAssignments();

        int splitPos = nextLevelToExpand(curAssignments);

        System.out.println("Splitting at Level: " + splitPos);
        System.out.println("*****************************");
        for (int rot : this.pruneMat.unprunedRCsAtPos(splitPos)) {
            int[] childAssignments = curAssignments.clone();
            childAssignments[splitPos] = rot;

            partFuncNode childNode = new partFuncNode(childAssignments);
            if (!nodeIsDefined(childNode)) {
                childNode.setLowerBoundLogZ(computeLowerBound(childNode));
//                childNode.setUpperBoundLogZ(computeUpperBound(childNode));
//                childNode.setScore(scoreNode(childNode));
                updateLowerBound(childNode);

                children.add(childNode);
            } else {//child node is leaf so we compute exact score
                double confE = getConfE(childNode);
                updateRunningSumLB(confE);
                numConfsEnumerated++;
            }
        }
        printEffectiveEpsilon();
        if (!isFullyAssigned(curNode)) {
            subtractUpperBound(node);
            //Now update upper bounds
            for (AStarNode childNode : children) {
                partFuncNode child = (partFuncNode) childNode;
                child.setUpperBoundLogZ(computeUpperBound(child));
                updateUpperBound(child);
                child.setScore(scoreNode(child));
            }
        }
        if (verbose) {
            System.out.println("LowerBound logZ: " + this.ef.log(this.lbZ.add(this.runningSum)).doubleValue());
            System.out.println("UpperBound logZ: " + this.ef.log(this.ubZ.add(this.runningSum)).doubleValue());
        }
        printEffectiveEpsilon();
        return children;
    }

    private int nextLevelToExpand(int[] partialConf) {
        if (useDynamicOrdering) {

            int bestLevel = -1;
            double bestLevelScore = Double.NEGATIVE_INFINITY;

            for (int level = 0; level < this.numPos; level++) {
                if (partialConf[level] < 0) {

                    double levelScore = scoreExpansionLevel(level, partialConf);

                    if (levelScore > bestLevelScore) {//higher score is better
                        bestLevelScore = levelScore;
                        bestLevel = level;
                    }
                }
            }

            if (bestLevel == -1) {
                throw new RuntimeException("ERROR: No next expansion level found for dynamic ordering");
            }

            return bestLevel;
        } else {//static ordering
            for (int level = 0; level < numPos; level++) {
                if (partialConf[level] < 0) {
                    return level;
                }
            }

            throw new RuntimeException("ERROR: Can't find next expansion level for fully defined conformation");
        }
    }

    private double scoreExpansionLevel(int level, int[] partialConf) {
        //We will score a level by new lower-bound
        //Thus the best level  is the level that most improves our lower bound

        int[] expandedConf = partialConf.clone();

        ArrayList<Integer> unprunedRCsAtLevel = this.pruneMat.unprunedRCsAtPos(level);
        ArrayList<Double> lowerBoundPerRC = new ArrayList<>();

        //Keep track of largest logZ lower bound for numerical accuracy (see below)
        double largestLowerBound = Double.NEGATIVE_INFINITY;

        for (int rc : unprunedRCsAtLevel) {
            expandedConf[level] = rc;
            double lowerBound = computeLowerBound(expandedConf);
            lowerBoundPerRC.add(lowerBound);
            largestLowerBound = Math.max(largestLowerBound, lowerBound);
        }

        //our new bound is sum e^bound for each RC
        //this is equivalent to e^largestBound * (sum e^(bound - largestBound)
        //which is helpful for numerical reasons (so we don't overflow)
        double score = 0.0;
        for (double bound : lowerBoundPerRC) {
            double normalizedBound = bound - largestLowerBound;
            score += Math.exp(normalizedBound);
        }
        score = Math.log(score) + largestLowerBound;

        return score;
    }

    private void subtractFromBounds(partFuncNode node) {
        this.lbZ = this.lbZ.subtract(this.ef.exp(node.getLowerBoundLogZ()));
        this.ubZ = this.ubZ.subtract(this.ef.exp(node.getUpperBoundLogZ()));
    }

    private void subtractLowerBound(partFuncNode node) {
        this.lbZ = this.lbZ.subtract(this.ef.exp(node.getLowerBoundLogZ()));
    }

    private void subtractUpperBound(partFuncNode node) {
        this.ubZ = this.ubZ.subtract(this.ef.exp(node.getUpperBoundLogZ()));
    }

    private void updateBounds(partFuncNode node) {
        this.lbZ = this.lbZ.add(this.ef.exp(node.getLowerBoundLogZ()));
        this.ubZ = this.ubZ.add(this.ef.exp(node.getUpperBoundLogZ()));
    }

    private void updateLowerBound(partFuncNode node) {
        this.lbZ = this.lbZ.add(this.ef.exp(node.getLowerBoundLogZ()));
    }

    private void updateUpperBound(partFuncNode node) {
        this.ubZ = this.ubZ.add(this.ef.exp(node.getUpperBoundLogZ()));
    }

    private void printEffectiveEpsilon() {
        double effectiveEpsilon = 1 - ((lbZ.add(this.runningSum)).divide((ubZ.add(this.runningSum)), this.ef.mc)).doubleValue();
        System.out.println("Effective Epsilon: " + effectiveEpsilon);
    }

    private void updateRunningSumLB(double confE) {
        //for decimal precision, we keep a normalizer, y, s.t. our sum is exp(y/constRT)*runningSum
        this.runningSum = this.runningSum.add(this.ef.exp(-confE / this.constRT));
    }

    private double getConfE(partFuncNode node) {
        if (!(nodeIsDefined(node))) {
            throw new RuntimeException("Cannot Compute energy for partial node");
        }
        //TODO: ADD functionality for EPIC/LUTE/ContinuousFlexibility
        return this.emat.getInternalEnergy(new RCTuple(node.getNodeAssignments()));
    }

    private boolean nodeIsDefined(partFuncNode node) {
        int[] nodeAssignments = node.getNodeAssignments();
        for (int rot : nodeAssignments) {
            if (rot == -1) {
                return false;
            }
        }
        return true;
    }

    @Override
    public boolean isFullyAssigned(AStarNode node) {
        double logLBZ = Math.log(this.lbZ.doubleValue() + this.runningSum.doubleValue());
        double logUBZ = Math.log(this.ubZ.doubleValue() + this.runningSum.doubleValue());
        return (logLBZ - logUBZ) >= Math.log(1 - this.epsilon);
    }

    @Override
    public AStarNode rootNode() {
        int[] conf = new int[this.numPos];
        Arrays.fill(conf, -1);//indicates the sequence is not assigned

        partFuncNode root = new partFuncNode(conf);
        root.setUpperBoundLogZ(computeUpperBound(root));
        root.setLowerBoundLogZ(computeLowerBound(root));
        root.setScore(scoreNode(root));
        updateBounds(root);
        printEffectiveEpsilon();
        return root;
    }

    public double computeEpsilonApprox(double epsilon) {
        this.epsilon = epsilon;
        this.nextConf();
        return Math.log(this.lbZ.doubleValue() + this.runningSum.doubleValue());
    }

    private double scoreNode(partFuncNode node) {
//        return -node.lbLogZ;
        return this.ef.exp(node.lbLogZ).subtract(this.ef.exp(node.ubLogZ)).doubleValue();
//        return mplp.optimizeMPLP(node.getNodeAssignments(), 1000);
    }

    private double computeLowerBound(partFuncNode node) {
        ReparamMRF mrf = new ReparamMRF(this.emat, this.pruneMat, node.getNodeAssignments(), this.eCut);
        SCMF_Clamp scmf = new SCMF_Clamp(mrf);
        double lbLogZ = scmf.getLogZLB();
        return lbLogZ;
    }

    private double computeLowerBound(int[] partialConf) {
        ReparamMRF mrf = new ReparamMRF(this.emat, this.pruneMat, partialConf, this.eCut);
        SCMF_Clamp scmf = new SCMF_Clamp(mrf);
        double lbLogZ = scmf.getLogZLB();
        return lbLogZ;
    }

    private double computeUpperBound(partFuncNode node) {
        if (useMapPert) {
            MapPerturbation mp = new MapPerturbation(this.emat, this.pruneMat);
            double ubLogZ = mp.calcUBLogZLPMax(node.getNodeAssignments(), this.numSamples);
            return ubLogZ;
        } else {
            /*            MarkovRandomField mrf = new MarkovRandomField(this.emat, this.pruneMat, node.getNodeAssignments(), this.eCut);
             TreeReweightedBeliefPropagation trbp = new TreeReweightedBeliefPropagation(mrf);
             double ubLogZ = trbp.getLogZ();
             */

//            double lbLogZ = computePartFunctionEstimate(emat, pruneMat, node.getNodeAssignments(), 10000);

            ReparamMRF rMRF = new ReparamMRF(this.emat, this.pruneMat, node.getNodeAssignments(), this.eCut);
            //TRBP2 trbp = new TRBP2(rMRF);
            TRBPSeq trbp = new TRBPSeq(rMRF);
            double ubLogZ2 = trbp.getLogZ();
            return ubLogZ2;
//            return getUpperBoundMPLP(node);
        }
    }
    
    private double getUpperBoundMPLP(partFuncNode node){
        double lpGMEC = mplp.optimizeMPLP(node.getNodeAssignments(), 1000);
        BigDecimal boltzmann = this.ef.exp(-lpGMEC/this.constRT);
        boltzmann = boltzmann.multiply(new BigDecimal(getNumConfsUnderNode(node.getNodeAssignments())));
        return this.ef.log(boltzmann).doubleValue();
    }
    private BigInteger getNumConfsUnderNode(int[] partialConf) {
        BigInteger numConfs = new BigInteger("1");
        for (int pos = 0; pos < partialConf.length; pos++) {
            if (partialConf[pos] == -1) {
                int numRCs = this.pruneMat.unprunedRCsAtPos(pos).size();
                numConfs = numConfs.multiply(new BigInteger(Integer.toString(numRCs)));
            }
        }
        return numConfs;
    }


    private double computePartFunctionEstimate(EnergyMatrix emat, PruningMatrix pruneMat, int[] partialConf, int numIter) {
        UpdatedPruningMatrix upm = new UpdatedPruningMatrix(pruneMat);
        for (int pos = 0; pos < partialConf.length; pos++) {
            if (partialConf[pos] != -1) {
                for (int rc : pruneMat.unprunedRCsAtPos(pos)) {
                    if (rc != partialConf[pos]) {
                        upm.markAsPruned(new RCTuple(pos, rc));
                    }
                }
            }
        }

        double gmecE = 0.0;
        double partFunc = 0.0;
        ConfTree tree = new ConfTree(emat, upm);
        for (int i = 0; i < numIter; i++) {
            int[] conf = tree.nextConf();
            if (conf == null) {
                break;
            }
            double energy = emat.getInternalEnergy(new RCTuple(conf));
            if (i == 0) {
                gmecE = -energy / constRT;
                partFunc += 1;
            } else {
                double normalized = (-energy / constRT) - gmecE;
                partFunc += Math.exp(normalized);
            }
        }
        return gmecE + Math.log(partFunc);

    }

}
