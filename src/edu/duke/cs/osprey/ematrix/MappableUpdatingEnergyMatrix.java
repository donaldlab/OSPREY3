package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.sharkstar.SHARKStarNodeScorer;

import java.util.List;

/**
 * Used to compute upper energy bound corrections in SHARK*
 *
 */
public class MappableUpdatingEnergyMatrix extends UpdatingEnergyMatrix<MappableTupE> {
    private SimpleConfSpace confSpace;
    private SHARKStarNodeScorer correctionUpperBoundHScorer;

    public MappableUpdatingEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target, ConfEnergyCalculator confECalc) {
        super(confSpace, target, confECalc);
        this.confSpace = confSpace;
        this.correctionUpperBoundHScorer = new SHARKStarNodeScorer(super.target, true);
    }

    public MappableUpdatingEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target) {
        super(confSpace, target);
        this.confSpace = confSpace;
        this.correctionUpperBoundHScorer = new SHARKStarNodeScorer(super.target, true);
    }

    @Override
    protected MappableTupE makeT(RCTuple tup, double val) {
        return null;
    }

    @Override
    protected TupETrie<MappableTupE> makeTrie(List<SimpleConfSpace.Position> positions) {
        return new MappableTupETrie(positions);
    }
}

