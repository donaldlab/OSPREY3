/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.partfunc;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.AStarTree;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.partitionfunctionbounds.MapPerturbation;
import edu.duke.cs.osprey.partitionfunctionbounds.MarkovRandomField;
import edu.duke.cs.osprey.partitionfunctionbounds.SelfConsistentMeanField;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ExpFunction;
import java.math.BigDecimal;
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

    BigDecimal runningSumLB;
    double runningSumExpNormalizer = 0.0;

    double epsilon = 0.1;

    int numPos;

    EnergyMatrix emat;
    PruningMatrix pruneMat;

    final double eCut = 0.0; //if we want to set an energy cutoff for MRFs

    boolean useMapPert = true;
    int numSamples = 10; //number of sample averages to take the max from;

    ExpFunction ef = new ExpFunction();
    double constRT = PoissonBoltzmannEnergy.constRT;

    boolean verbose = true;
    
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

    }

    private void init(SearchProblem sp, PruningMatrix aPruneMat, boolean useEPIC) {
        numPos = sp.confSpace.numPos;
        this.lbZ = new BigDecimal("0.0");
        this.ubZ = new BigDecimal("0.0");
        this.runningSumLB = new BigDecimal("0.0");

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
    }

    @Override
    public ArrayList<AStarNode> getChildren(AStarNode curNode) {
        partFuncNode node = (partFuncNode) curNode;
        ArrayList<AStarNode> children = new ArrayList<>();

        subtractFromBounds(node);
        
        int[] curAssignments = node.getNodeAssignments();
        for (int splitPos = 0; splitPos < this.numPos; splitPos++) {
            if (curAssignments[splitPos] < 0) {
                for (int rot : this.pruneMat.unprunedRCsAtPos(splitPos)) {
                    int[] childAssignments = curAssignments.clone();
                    childAssignments[splitPos] = rot;

                    partFuncNode childNode = new partFuncNode(childAssignments);
                    if (!nodeIsDefined(childNode)) {
                        childNode.setLowerBoundLogZ(computeLowerBound(childNode));
                        childNode.setUpperBoundLogZ(computeUpperBound(childNode));
                        childNode.setScore(scoreNode(childNode));
                        updateBounds(childNode);

                        children.add(childNode);
                    } else {//child node is leaf so we compute exact score
                        double confE = getConfE(childNode);
                        updateRunningSumLB(confE);
                    }
                }
                return children;
            }
        }
        throw new RuntimeException("ERROR: Not splittable position found");
    }

    private void subtractFromBounds(partFuncNode node){
        this.lbZ = this.lbZ.subtract(this.ef.exp(node.getLowerBoundLogZ()));
        this.ubZ = this.ubZ.subtract(this.ef.exp(node.getUpperBoundLogZ()));
    }
    
    private void updateBounds(partFuncNode node){
        this.lbZ = this.lbZ.add(this.ef.exp(node.getLowerBoundLogZ()));
        this.ubZ = this.ubZ.add(this.ef.exp(node.getUpperBoundLogZ()));
        if (verbose){
            System.out.println("LowerBound Z: "+(this.lbZ.doubleValue() + this.runningSumLB.doubleValue()));
            System.out.println("UpperBound Z: "+this.ubZ.doubleValue());
        }
        double effectiveEpsilon = 1-(lbZ.doubleValue()/(ubZ.doubleValue()));
        System.out.println("Effective Epsilon: "+effectiveEpsilon);
    }
    
    private void updateRunningSumLB(double confE) {
        //for decimal precision, we keep a normalizer, y, s.t. our sum is exp(y/constRT)*runningSum
        this.runningSumLB = this.runningSumLB.add(this.ef.exp(-confE/this.constRT));
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
        double logLBZ = Math.log(this.lbZ.doubleValue() + this.runningSumLB.doubleValue());
        double logUBZ = Math.log(this.ubZ.doubleValue());
        return (logLBZ-logUBZ) >= Math.log(1-this.epsilon);
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
        
        return root;
    }

    public double computeEpsilonApprox(double epsilon){
        this.epsilon = epsilon;
        this.nextConf();
        return this.lbZ.doubleValue()+this.runningSumLB.doubleValue();
    }
    
    private double scoreNode(partFuncNode node) {
        return node.getLowerBoundLogZ() - node.getUpperBoundLogZ();
    }

    private double computeLowerBound(partFuncNode node) {
        MarkovRandomField mrf = new MarkovRandomField(this.emat, this.pruneMat, node.getNodeAssignments(), this.eCut);
        SelfConsistentMeanField scmf = new SelfConsistentMeanField(mrf);
        scmf.run();
        double lbLogZ = scmf.calcLBLogZ();
        return lbLogZ;
    }

    private double computeUpperBound(partFuncNode node) {
        if (useMapPert) {
            MapPerturbation mp = new MapPerturbation(this.emat, this.pruneMat);
            double ubLogZ = mp.calcUBLogZLPMax(node.getNodeAssignments(), this.numSamples);
            return ubLogZ;
        } else {
            throw new RuntimeException("ERROR: PartFuncTree -- MapPert is only option for upper bound");
        }
    }

}
