package edu.duke.cs.osprey.energy.forcefield.amber;

import java.util.List;

public record BondLengthParameter(AtomSymbolAndMass IBT, AtomSymbolAndMass JBT, float RK, float REQ) implements HasAtoms {

    @Override
    public List<AtomSymbolAndMass> atoms() {
        return List.of(IBT, JBT);
    }

}
