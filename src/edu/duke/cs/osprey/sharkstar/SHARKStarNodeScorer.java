package edu.duke.cs.osprey.sharkstar;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;

import static edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarBound.confMatch;
import static edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarBound.debugConf;

public class SHARKStarNodeScorer implements AStarScorer {
    protected EnergyMatrix emat;
    protected MathTools.Optimizer opt = MathTools.Optimizer.Minimize;
    /*
    protected int[] debugConf;
     */

    protected int[][][] bestPairs;

    public SHARKStarNodeScorer(EnergyMatrix emat, boolean negated) {
        this.emat = emat;
        if(negated)
            opt = MathTools.Optimizer.Maximize;
        /*
        this.debugConf = new int[]{};

         */
        init();
    }

    protected void init(){
        int maxNumRots = Arrays.stream(emat.getNumConfAtPos())
                .boxed()
                .max(Integer::compare)
                .get();
        bestPairs = new int[emat.getNumPos()][maxNumRots][emat.getNumPos()];
        fillBestPairMatrix();
    }

    private void fillBestPairMatrix(){
        //TODO: have this consider multiple sequences?
        for (int pos = 0; pos < emat.getNumPos(); pos++) {
            for (int rc = 0; rc < emat.getNumConfAtPos(pos); rc++){
                for (int partnerPos = 0; partnerPos < emat.getNumPos(); partnerPos++){
                    if(pos == partnerPos)
                        continue;
                    int bestPartnerRC = 0;
                    double bestPairwiseEnergy = opt.initDouble();
                    for(int partnerRC = 0; partnerRC < emat.getNumConfAtPos(partnerPos); partnerRC++){
                       double partnerEnergy = emat.getEnergy(pos, rc, partnerPos, partnerRC);
                       if( opt.isBetter(partnerEnergy, bestPairwiseEnergy)){
                           bestPartnerRC = partnerRC;
                           bestPairwiseEnergy = partnerEnergy;
                       }
                   }
                   bestPairs[pos][rc][partnerPos] = bestPartnerRC;
                }
            }
        }
    }

    @Override
    public AStarScorer make() {
        return new SHARKStarNodeScorer(emat, false);
    }

    public double calc(ConfIndex confIndex, Sequence seq, SimpleConfSpace confSpace) {
        return calc(confIndex, seq.makeRCs(confSpace));
    }

    public BigDecimal calcGScore(ConfIndex confIndex, RCs rcs){
        /** calcGScore
         *
         * Calculates the contribution to the partition function bound of the defined residues ONLY.
         *
         * Intended primarily for debugging upperBound error - GTH 200114
         */

        BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
        BigDecimal pfuncBound = BigDecimal.ONE;

        // defined one-body energies
        for (int i=0; i<confIndex.numDefined; i++) {
            int pos1 = confIndex.definedPos[i];
            int rc1 = confIndex.definedRCs[i];
            pfuncBound = pfuncBound.multiply(bcalc.calc(emat.getOneBody(pos1, rc1)));
        }

        // defined pairwise energies
        for (int i=0; i<confIndex.numDefined; i++) {
            int pos1 = confIndex.definedPos[i];
            int rc1 = confIndex.definedRCs[i];

            for (int j=0; j<i; j++) {
                int pos2 = confIndex.definedPos[j];
                int rc2 = confIndex.definedRCs[j];

                pfuncBound = pfuncBound.multiply(bcalc.calc(emat.getPairwise(pos1, rc1, pos2, rc2)));
            }
        }
        return pfuncBound;
    }

    public BigDecimal calcCombinedScore(ConfIndex confIndex, RCs rcs){
        /**calcCombinedScore
         *
         * Calculates the bound on the partition function value at a node,
         * bounding the sum of boltzmann-weighted energies of all conformations
         * that are compatible with the given conformation.
         *
         * Intended primarily for debugging upperBound error - GTH 200114
         */

        BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
        BigDecimal pfuncBound = BigDecimal.ONE;

        // defined one-body energies
        for (int i=0; i<confIndex.numDefined; i++) {
            int pos1 = confIndex.definedPos[i];
            int rc1 = confIndex.definedRCs[i];
            pfuncBound = pfuncBound.multiply(bcalc.calc(emat.getOneBody(pos1, rc1)));
        }

        // defined pairwise energies
        for (int i=0; i<confIndex.numDefined; i++) {
            int pos1 = confIndex.definedPos[i];
            int rc1 = confIndex.definedRCs[i];

            for (int j=0; j<i; j++) {
                int pos2 = confIndex.definedPos[j];
                int rc2 = confIndex.definedRCs[j];

                pfuncBound = pfuncBound.multiply(bcalc.calc(emat.getPairwise(pos1, rc1, pos2, rc2)));
            }
        }

        // now we iterate through the undefined residues
        for (int undefinedPosIndex1 = 0; undefinedPosIndex1 < confIndex.numUndefined; undefinedPosIndex1++) {
            int undefinedPos1 = confIndex.undefinedPos[undefinedPosIndex1];
            BigDecimal residueSum = BigDecimal.ZERO;
            for (int rot1 : rcs.get(undefinedPos1)) {
                double rotEnergy = emat.getEnergy(undefinedPos1, rot1); // get the singletons

                for (int definedPosIndex = 0; definedPosIndex < confIndex.numDefined; definedPosIndex ++) { // iterate through the defined residues
                    int definedPos = confIndex.definedPos[definedPosIndex];
                    int definedRC = confIndex.definedRCs[definedPosIndex];
                    rotEnergy += emat.getEnergy(undefinedPos1, rot1, definedPos, definedRC); // get the defined pairwise
                }
                for (int undefinedPosIndex2 = 0; undefinedPosIndex2 < confIndex.numUndefined; undefinedPosIndex2++) { // iterate through the undefined residues
                    int undefinedPos2 = confIndex.undefinedPos[undefinedPosIndex2];
                    if (undefinedPos2 >= undefinedPos1)
                        continue;
                    rotEnergy+=emat.getEnergy(undefinedPos1, rot1, undefinedPos2, bestPairs[undefinedPos1][rot1][undefinedPos2]);
                }
                residueSum = residueSum.add(bcalc.calc(rotEnergy), bcalc.mathContext);
            }
            pfuncBound = pfuncBound.multiply(residueSum, PartitionFunction.decimalPrecision);
        }

        return pfuncBound;
    }

    /* Assumes: that the rcs contain only the sequence in question. In this case, we need only
     *  sum over all unassigned positions. Returns a lower bound on the ensemble energy.
     *  Note: I currently exponentiate and log for compatibilty. This could be optimized.*/
    @Override
    public double calc(ConfIndex confIndex, RCs rcs) {
        BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
        double freeEnergyBound = 0.0;
        BigDecimal pfuncBound = BigDecimal.ONE;

        int[] conf = confIndex.makeConf();
        if(confMatch(conf, debugConf)){
            System.out.println("Node Conf: "+conf.toString());
            RCTuple testTuple = new RCTuple(conf);
            double energyCheck = emat.getInternalEnergy(testTuple);
            System.out.println("Energy of " + testTuple + ": " + energyCheck);
            System.out.println("Gotcha-calc0");
        }
        for (int undefinedPosIndex1 = 0; undefinedPosIndex1 < confIndex.numUndefined; undefinedPosIndex1++) {
            int undefinedPos1 = confIndex.undefinedPos[undefinedPosIndex1];
            BigDecimal residueSum = BigDecimal.ZERO;
            // ArrayList to use the log-sum-exp trick
            ArrayList <Double> rotEnergies = new ArrayList();

            for (int rot1 : rcs.get(undefinedPos1)) {
                //Iterate over the residues, computing bounds on the log partition function contributions...
                double rotEnergy = emat.getEnergy(undefinedPos1, rot1); // get the singletons

                for (int definedPosIndex = 0; definedPosIndex < confIndex.numDefined; definedPosIndex ++) { // iterate through the defined residues
                    int definedPos = confIndex.definedPos[definedPosIndex];
                    int definedRC = confIndex.definedRCs[definedPosIndex];
                    rotEnergy += emat.getEnergy(undefinedPos1, rot1, definedPos, definedRC); // get the defined pairwise
                }
                for (int undefinedPosIndex2 = 0; undefinedPosIndex2 < confIndex.numUndefined; undefinedPosIndex2++) { // iterate through the undefined residues
                    int undefinedPos2 = confIndex.undefinedPos[undefinedPosIndex2];
                    if (undefinedPos2 >= undefinedPos1)
                        continue;
                    rotEnergy+=emat.getEnergy(undefinedPos1, rot1, undefinedPos2, bestPairs[undefinedPos1][rot1][undefinedPos2]);
                }
                if(confMatch(conf, debugConf)){
                    System.out.println("Gotcha-calc");
                    int[] conf2 = confIndex.makeConf();
                    conf2[undefinedPos1] = rot1;
                    RCTuple testTuple = new RCTuple(conf2);
                    double energyCheck = emat.getInternalEnergy(testTuple);
                    System.out.println("Child Conf: "+conf2.toString());
                    System.out.println("Energy of " + testTuple + ": " + energyCheck);
                    System.out.println("Rot energy of "+ testTuple + ": " + rotEnergy);
                }
                // Store the rotamer energy for later residue energy calculation
                rotEnergies.add(rotEnergy);
                //OLD
                residueSum = residueSum.add(bcalc.calc(rotEnergy),bcalc.mathContext);
            }
            // Compute the residue sum via logSumExp
            double residueFreeEnergy = bcalc.logSumExp(rotEnergies);
            if(confMatch(conf, debugConf)){
                System.out.println("Gotcha-calc2");
                System.out.println("End residue sum: "+residueSum);
                System.out.println("End residue sum (via logsumexp): "+residueFreeEnergy);
                System.out.println("bcalc residue sum (via logsumexp): "+bcalc.calc(residueFreeEnergy));
            }
            //OLD
            pfuncBound = pfuncBound.multiply(residueSum, PartitionFunction.decimalPrecision);
            // Add the residue sum to total bound
            freeEnergyBound += residueFreeEnergy;
        }
        if(confMatch(conf, debugConf)){
            System.out.println("Gotcha-calc3");
            System.out.printf("End bound: %12.4e\n", pfuncBound);
            System.out.println("End \'energy\': "+bcalc.freeEnergy(pfuncBound));
            System.out.println("End \'energy\' (via logsumexp): "+freeEnergyBound);
        }
        //return bcalc.freeEnergy(pfuncBound);
        return freeEnergyBound;
    }

    /*
    public void setDebugConf(int[] conf){
        this.debugConf = conf;
    }

     */
}
