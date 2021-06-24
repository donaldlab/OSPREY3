package edu.duke.cs.osprey.energy.forcefield.amber;

import java.util.List;

public record ImproperDihederalParameter(AtomSymbolAndMass IPT, AtomSymbolAndMass JPT, AtomSymbolAndMass KPT,
                                         AtomSymbolAndMass LPT, float PK, float PHASE, float PN) implements HasAtoms {
    @Override
    public List<AtomSymbolAndMass> atoms() {
        return List.of(IPT, JPT, KPT, LPT);
    }
}
