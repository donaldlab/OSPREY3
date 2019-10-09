package edu.duke.cs.osprey.sharkstar;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;

import static edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarBound.confMatch;
import static edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarBound.debugConf;

public class SHARKStarNodeScorer implements AStarScorer {
    protected EnergyMatrix emat;
    protected MathTools.Optimizer opt = MathTools.Optimizer.Minimize;

    public SHARKStarNodeScorer(EnergyMatrix emat, boolean negated) {
        this.emat = emat;
        if(negated)
            opt = MathTools.Optimizer.Maximize;
    }

    @Override
    public AStarScorer make() {
        return new SHARKStarNodeScorer(emat, false);
    }

    public double calc(ConfIndex confIndex, Sequence seq, SimpleConfSpace confSpace) {
        return calc(confIndex, seq.makeRCs(confSpace));
    }

    /* Assumes: that the rcs contain only the sequence in question. In this case, we need only
     *  sum over all unassigned positions. Returns a lower bound on the ensemble energy.
     *  Note: I currently exponentiate and log for compatibilty. This could be optimized.*/
    @Override
    public double calc(ConfIndex confIndex, RCs rcs) {
        BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
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
            for (int rot1 : rcs.get(undefinedPos1)) {
                double rotEnergy = emat.getEnergy(undefinedPos1, rot1);
                for (int definedPosIndex = 0; definedPosIndex < confIndex.numDefined; definedPosIndex ++) {
                    int definedPos = confIndex.definedPos[definedPosIndex];
                    int definedRC = confIndex.definedRCs[definedPosIndex];
                    rotEnergy += emat.getEnergy(undefinedPos1, rot1, definedPos, definedRC);
                }
                for (int undefinedPosIndex2 = 0; undefinedPosIndex2 < confIndex.numUndefined; undefinedPosIndex2++) {
                    int undefinedPos2 = confIndex.undefinedPos[undefinedPosIndex2];
                    if (undefinedPos2 >= undefinedPos1)
                        continue;
                    double bestPair = Double.MAX_VALUE;
                    for (int rot2 : rcs.get(undefinedPos2)) {
                        bestPair = opt.opt(bestPair, emat.getEnergy(undefinedPos1, rot1, undefinedPos2, rot2));
                    }
                    rotEnergy+= bestPair;
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
                residueSum = residueSum.add(bcalc.calc(rotEnergy));
            }
            if(confMatch(conf, debugConf)){
                System.out.println("Gotcha-calc2");
                System.out.println("End residue sum: "+residueSum);
            }
            pfuncBound = pfuncBound.multiply(residueSum, PartitionFunction.decimalPrecision);
        }
        if(confMatch(conf, debugConf)){
            System.out.println("Gotcha-calc3");
            System.out.println("End bound: "+bcalc.freeEnergy(pfuncBound));
        }
        return bcalc.freeEnergy(pfuncBound);
    }
}
