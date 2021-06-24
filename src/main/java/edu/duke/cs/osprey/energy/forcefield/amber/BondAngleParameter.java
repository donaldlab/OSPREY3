package edu.duke.cs.osprey.energy.forcefield.amber;

import java.util.List;

public record BondAngleParameter(AtomSymbolAndMass ITT, AtomSymbolAndMass JTT, AtomSymbolAndMass KTT, float TK, float TEQ) implements HasAtoms {
    @Override
    public List<AtomSymbolAndMass> atoms() {
        return List.of(ITT, JTT, KTT);
    }
}
