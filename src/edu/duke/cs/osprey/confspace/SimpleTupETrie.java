package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TupE;
import edu.duke.cs.osprey.confspace.TupETrie;

import java.util.List;

public class SimpleTupETrie extends TupETrie<TupE> {

    public SimpleTupETrie(List<SimpleConfSpace.Position> positions) {
        super(positions);
    }

    @Override
    protected TupE makeT(String repr) {
        return new TupE(repr);
    }
}
