package edu.duke.cs.osprey.energy.forcefield.amber;

import java.util.List;

public record Six12PotentialCoefficient(AtomSymbolAndMass LTYNB, float A, float C) implements HasAtoms {
    @Override
    public List<AtomSymbolAndMass> atoms() {
        return List.of(LTYNB);
    }
}
