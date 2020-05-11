package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;

import java.util.List;

public class SimpleUpdatingEnergyMatrix extends UpdatingEnergyMatrix<TupE> {
    public SimpleUpdatingEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target, ConfEnergyCalculator confECalc) {
        super(confSpace, target, confECalc);
    }

    public SimpleUpdatingEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target) {
        super(confSpace, target);
    }

    @Override
    protected TupE makeT(RCTuple tup, double val) {
        return new TupE(tup, val);
    }

    @Override
    protected TupETrie<TupE> makeTrie(List<SimpleConfSpace.Position> positions) {
        return new SimpleTupETrie(positions);
    }
}
