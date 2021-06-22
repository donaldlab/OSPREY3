package edu.duke.cs.osprey.energy.forcefield.amber;

import java.util.List;

public record HBond10_12PotentialParameter(AtomSymbolAndMass KT1, AtomSymbolAndMass KT2, float A,
                                           float B/*, float ASOLN, float BSOLN, float HCUT, int IC*/) implements HasAtoms {
    @Override
    public List<AtomSymbolAndMass> atoms() {
        return List.of(KT1, KT2);
    }
} // These other fields, while present in standard, are not used
