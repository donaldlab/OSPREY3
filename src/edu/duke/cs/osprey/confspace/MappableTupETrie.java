package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.confspace.MappableTupE;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TupETrie;

import java.util.List;

public class MappableTupETrie extends TupETrie<MappableTupE> {
    public MappableTupETrie(List<SimpleConfSpace.Position> positions) {
        super(positions);
    }

    @Override
    protected MappableTupE makeT(String repr) {
        return new MappableTupE(repr);
    }
}
