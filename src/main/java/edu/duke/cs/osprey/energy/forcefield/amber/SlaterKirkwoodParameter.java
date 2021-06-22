package edu.duke.cs.osprey.energy.forcefield.amber;

import java.util.List;

public record SlaterKirkwoodParameter(AtomSymbolAndMass LTYNB, float POL, float XNEFF, float RMIN) implements HasAtoms {
    @Override
    public List<AtomSymbolAndMass> atoms() {
        return List.of(LTYNB);
    }
}
