package edu.duke.cs.osprey.energy.forcefield.amber;

import java.util.List;

public record VanDerWaalsRadius(AtomSymbolAndMass LTYNB, float R, float EDEP) implements HasAtoms {
    @Override
    public List<AtomSymbolAndMass> atoms() {
        return List.of(LTYNB);
    }
}
