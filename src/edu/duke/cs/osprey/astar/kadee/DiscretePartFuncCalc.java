/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.kadee;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ExpFunction;
import java.math.BigDecimal;
import java.math.BigInteger;

/**
 * This class is a quick way to quickly estimate a discrete partition function
 * using an A* tree. It is not meant for production. It is just here for testing
 * purposes.
 *
 * @author hunternisonoff
 */
public class DiscretePartFuncCalc {

    BigDecimal Z = new BigDecimal(0.0);

    EnergyMatrix emat;
    PruningMatrix pruneMat;

    ConfTree confTree;

    int numConfsToEnumerate;
    double epsilon; 
    
    boolean partFuncIsCalculated = false;
    
    final ExpFunction ef = new ExpFunction();
    final double constRT = PoissonBoltzmannEnergy.constRT;

    boolean verbose = false;
    
    public DiscretePartFuncCalc(EnergyMatrix emat, PruningMatrix pruneMat, int numConfs, double epsilon) {
        //Check to make sure we don't have any null objects
        if (emat == null) {
            throw new NullPointerException("Energy Matrix is Null");
        }
        if (pruneMat == null) {
            throw new NullPointerException("Pruning Matrix is Null");
        }
        if (numConfs <= 0){
            throw new RuntimeException("Error: Number of Confs to Enumerate must be greater than 0");
        }
        if (epsilon < 0 || epsilon > 1){
            throw new RuntimeException("Error: epsilon must be between 0 and 1. It is currently: "+epsilon);
        }
        //All clear: Initialize our Conf Tree
        this.emat = emat;
        this.pruneMat = pruneMat;
        this.confTree = new ConfTree(emat, pruneMat);
        this.epsilon = epsilon;
    }

    /**
     * Computes the partition function
     */
    public void calcPartFunc() {
        //First we get the gmec
        //This will help us with numerical accuracy by factoring out the Boltzmann-
        //weighted gmec energy in our sum
        int[] gmec = confTree.nextConf();
        double gmecE = emat.getInternalEnergy(new RCTuple(gmec));
        this.Z = this.Z.add(new BigDecimal(BigInteger.ONE)); //Z = e^(-gmecE/RT)(1 + ...)
        int confsEnumurated = 1; //(the gmec)
        boolean converged = false;
        while (!converged) {
            int[] conf = confTree.nextConf();
            double confE = emat.getInternalEnergy(new RCTuple(conf));
            double deltaE = gmecE - confE;
            this.Z = this.Z.add(this.ef.exp(-deltaE / this.constRT));
            confsEnumurated++;
            //check convergence criteria
            converged = hasConverged(deltaE, confsEnumurated);
        }
        partFuncIsCalculated = true;
    }

    /**
     * Computes the natural log of the partition function
     * @return 
     */
    public double calcLogZ(){
        if (!partFuncIsCalculated){
            calcPartFunc();
        }
        return ef.logToDouble(Z);
    }
    
    /**
     * Checks if we are finished calculating the partition function
     * @param currentConfE the energy of the last conformation to come off the tree
     * @param numConfsEnumerated the number of conformations that we have enumerated so far
     * @return 
     */
    private boolean hasConverged(double currentConfE, int numConfsEnumerated) {
        BigDecimal lowerBound = this.Z;
        BigDecimal upperBound = this.ef.exp(-currentConfE/this.constRT).multiply(new BigDecimal(calcNumRemainingConfs()));
        if (verbose){
            System.out.println("Lower Bound Z: "+lowerBound.doubleValue());
            System.out.println("Upper Bound Z: "+ upperBound.doubleValue());
            System.out.println("Gap between Bounds: "+(upperBound.doubleValue() - lowerBound.doubleValue()));
        }
        if (upperBound.doubleValue() - lowerBound.doubleValue() < epsilon) {
            return true;
        }
        else if (numConfsEnumerated >= this.numConfsToEnumerate){
            return true;
        }
        return false;
    }

    /**
     * Calculated the total number of remaining confs (or a rough estimate)
     * This is done by adding the size of the search problems of all the nodes in the queue
     * @return 
     */
    private BigInteger calcNumRemainingConfs(){
        BigInteger numRemainingConfs = new BigInteger("0");
        for (AStarNode node : this.confTree.pq){
            numRemainingConfs = numRemainingConfs.add(calcSearchSpace(node));
        }
        return numRemainingConfs;
    }

    /**
     * Calculates the size of the search space underneath a particular node in the tree;
     * @param node the node 
     * @return 
     */
    private BigInteger calcSearchSpace(AStarNode node) {
        BigInteger searchSpaceSize = new BigInteger("1");
        int[] assignments = node.getNodeAssignments();
        for (int pos : assignments){
            if (assignments[pos] == -1){
                int numUnprunedAtPos = pruneMat.unprunedRCsAtPos(pos).size();
                searchSpaceSize = searchSpaceSize.multiply(new BigInteger(Integer.toString(numUnprunedAtPos)));
            }
        }
        return searchSpaceSize;
    }
    
}
