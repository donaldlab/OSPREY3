package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class ManualOrder implements AStarOrder {

    private AStarScorer gscorer;
    private AStarScorer hscorer;

    private List<Integer> posOrder;

    public ManualOrder(){}

    @Override
    public void setScorers(AStarScorer gscorer, AStarScorer hscorer) {
        this.gscorer = gscorer;
        this.hscorer = hscorer;
    }

    @Override
    public boolean isDynamic() {
        return false;
    }

    @Override
    public int getNextPos(ConfIndex confIndex, RCs rcs) {
        if (posOrder == null) {
            posOrder = calcPosOrder(confIndex, rcs);
        }
        return posOrder.get(confIndex.node.getLevel());
    }

    private List<Integer> calcPosOrder(ConfIndex confIndex, RCs rcs){
        // init permutation array with only undefined positions and score them
        List<Integer> undefinedOrder = new ArrayList<>();
        for (int posi=0; posi<confIndex.numUndefined; posi++) {
            int pos = confIndex.undefinedPos[posi];
            undefinedOrder.add(pos);
        }
        // Do not sort, as we want the residues to remain in the order they are in
        // prepend the defined positions to build the full order
        List<Integer> order = new ArrayList<>();
        for (int posi=0; posi<confIndex.numDefined; posi++) {
            int pos = confIndex.definedPos[posi];
            order.add(pos);
        }
        order.addAll(undefinedOrder);

        return order;
    }

}
