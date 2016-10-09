/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ExpFunction;

import java.math.BigDecimal;
import java.math.BigInteger;

/**
 * @author hmn5
 */
public class DiscretePartFunc {

    double epsilon;
    EnergyMatrix emat;
    PruningMatrix pruneMat;

    double constRT = PoissonBoltzmannEnergy.constRT;
    ExpFunction ef = new ExpFunction();

    boolean printEffectiveEpsilon = true;

    double logZ;

    public long totalTime;
    public double effectiveEpsilonReached;
    public boolean finishedInTime;
    double lowerBoundLogZ;
    double upperBoundLogZ;

    public DiscretePartFunc(SearchProblem sp, double epsilon) {
        if (epsilon < 0 || epsilon >= 1.0) {
            throw new RuntimeException("Epsilon must be between 0 and 1");
        }
        this.epsilon = epsilon;
        this.emat = sp.emat;
        this.pruneMat = sp.pruneMat;

        this.logZ = computePartFunc();
    }

    public DiscretePartFunc(EnergyMatrix emat, PruningMatrix pruneMat, double epsilon) {
        this.epsilon = epsilon;
        this.emat = emat;
        this.pruneMat = pruneMat;

        this.logZ = computePartFunc();
    }

    public DiscretePartFunc(EnergyMatrix emat, PruningMatrix pruneMat, double epsilon, double maxTimeLimitMilli) {
        this.epsilon = epsilon;
        this.emat = emat;
        this.pruneMat = pruneMat;

        this.logZ = computePartFunc(maxTimeLimitMilli);
    }

    public double getLogZ() {
        return this.logZ;
    }

    double computePartFunc() {
        ConfTree tree = new ConfTree(emat, pruneMat);
        double normalize;
        double runningSum = 1.0;

        BigInteger numConfsRemaining = getNumConfs(pruneMat);

        int[] gmecConf = tree.nextConf();
        numConfsRemaining = numConfsRemaining.subtract(new BigInteger("1"));

        normalize = (-emat.getInternalEnergy(new RCTuple(gmecConf))) / constRT;
        boolean epsilonReached = false;
        boolean epsilonCannotBeReached = false;
        while (!epsilonReached && !epsilonCannotBeReached) {
            int[] nextConf = tree.nextConf();
            numConfsRemaining = numConfsRemaining.subtract(new BigInteger("1"));
            if (nextConf == null) {
                System.out.println("Conf was null");
                epsilonCannotBeReached = true;
            } else {
                double confE = (-emat.getInternalEnergy(new RCTuple(nextConf))) / constRT;
                double normalizedE = (confE - normalize);
                runningSum += Math.exp(normalizedE);
                BigDecimal upperBound = getUpperBound(runningSum, normalizedE, numConfsRemaining);
                BigDecimal lowerBound = new BigDecimal(runningSum);

                double effectiveEpsilon = getEffectiveEpsilon(upperBound, lowerBound);
                if (printEffectiveEpsilon) {
                    System.out.println("Effective Epsilon: " + effectiveEpsilon + "  logZ: " + (normalize + Math.log(runningSum)));
                }
                if (effectiveEpsilon <= epsilon) {
                    epsilonReached = true;
                }
            }
        }
        return (normalize) + Math.log(runningSum);
    }

    double computePartFunc(double maxTimeMilli) {
        long startTime = System.currentTimeMillis();
        ConfTree tree = new ConfTree(emat, pruneMat);
        double normalize;
        double runningSum = 1.0;

        BigInteger numConfsRemaining = getNumConfs(pruneMat);

        int[] gmecConf = tree.nextConf();
        numConfsRemaining = numConfsRemaining.subtract(new BigInteger("1"));

        normalize = (-emat.getInternalEnergy(new RCTuple(gmecConf))) / constRT;
        boolean epsilonReached = false;
        boolean epsilonCannotBeReached = false;

        double numConfs = 0;
        while (!epsilonReached && !epsilonCannotBeReached) {
            int[] nextConf = tree.nextConf();
            numConfsRemaining = numConfsRemaining.subtract(new BigInteger("1"));
            if (nextConf == null) {
                System.out.println("Conf was null");
                epsilonCannotBeReached = true;
            } else {
                double confE = (-emat.getInternalEnergy(new RCTuple(nextConf))) / constRT;
                double normalizedE = (confE - normalize);
                runningSum += Math.exp(normalizedE);
                BigDecimal upperBound = getUpperBound(runningSum, normalizedE, numConfsRemaining);
                BigDecimal lowerBound = new BigDecimal(runningSum);

                if (numConfs % 100 == 0) {
                    double effectiveEpsilon = getEffectiveEpsilon(upperBound, lowerBound);
                    if (printEffectiveEpsilon) {
                        System.out.println("Effective Epsilon: " + effectiveEpsilon + "  logZ: " + (normalize + Math.log(runningSum)));
                    }
                    if (effectiveEpsilon <= epsilon) {
                        epsilonReached = true;
                        this.lowerBoundLogZ = (normalize) + ef.log(lowerBound).doubleValue();
                        this.upperBoundLogZ = (normalize) + ef.log(upperBound).doubleValue();
                    } else if (System.currentTimeMillis() - startTime > maxTimeMilli) {
                        System.out.println("Effective Epsilon after " + maxTimeMilli + " milliseconds: " + effectiveEpsilon);
                        this.effectiveEpsilonReached = effectiveEpsilon;
                        finishedInTime = false;
                        this.lowerBoundLogZ = (normalize) + ef.log(lowerBound).doubleValue();
                        this.upperBoundLogZ = (normalize) + ef.log(upperBound).doubleValue();
                        return (normalize) + Math.log(runningSum);
                    }
                }
            }
        }
        this.totalTime = System.currentTimeMillis() - startTime;
        finishedInTime = true;
        return (normalize) + Math.log(runningSum);
    }

    public double getLowerBoundLogZ() {
        return lowerBoundLogZ;
    }

    public double getUpperBoundLogZ() {
        return upperBoundLogZ;
    }

    double getEffectiveEpsilon(BigDecimal upperBound, BigDecimal lowerBound) {
        return 1 - (lowerBound.divide(upperBound, ef.mc).doubleValue());
    }

    BigDecimal getUpperBound(double runningSum, double lastNormalizeE, BigInteger numConfsRemaining) {
        BigDecimal lastConfBW = new BigDecimal(Math.exp(lastNormalizeE));
        BigDecimal numConfsRemainingD = new BigDecimal(numConfsRemaining);
        BigDecimal boundOnRemaining = lastConfBW.multiply(numConfsRemainingD);
        BigDecimal runningSumD = new BigDecimal(runningSum);

        BigDecimal upperBound = runningSumD.add(boundOnRemaining);
        return upperBound;
    }

    BigInteger getNumConfs(PruningMatrix pruneMat) {
        BigInteger numConfs = new BigInteger("1");
        for (int pos = 0; pos < pruneMat.numPos(); pos++) {
            int numUnpruned = pruneMat.unprunedRCsAtPos(pos).size();
            numConfs = numConfs.multiply(new BigInteger(Integer.toString(numUnpruned)));
        }
        return numConfs;
    }

}
