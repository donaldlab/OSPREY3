package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.util.*;
import java.util.stream.Collectors;

/** RespectfulSBLDOrder
 *
 * Implements the standard StaticBiggestLowerboundDifferenceOrder scheme, with the caveat that
 * residues from different strands cannot be interleaved.
 *
 * Strands will be ordered by the order in which they were given to the confspace.
 *
 * This is primarily useful for energy landscape visualization with MARK*.
 */
public class RespectfulSBLDOrder extends StaticBiggestLowerboundDifferenceOrder{

    private List<Integer> posOrder;
    private SimpleConfSpace confSpace;

    public RespectfulSBLDOrder(SimpleConfSpace confSpace){
        this.confSpace = confSpace;
    }

    @Override
    public int getNextPos(ConfIndex confIndex, RCs rcs) {
        // if there is only one strand, default to SBLDOrder
        if (confSpace.strands.size() <= 1){
            posOrder = super.calcPosOrder(confIndex, rcs);
        }else{
            posOrder = calcPosOrder(confIndex, rcs, confSpace);
        }
        return posOrder.get(confIndex.node.getLevel());
    }
    private List<Integer> calcPosOrder(ConfIndex confIndex, RCs rcs, SimpleConfSpace confspace) {
        // init permutation array with only undefined positions and score them
        List<Integer> undefinedOrder = new ArrayList<Integer>();
        Map<Integer, Double > scores = new TreeMap<>();
        for (int posi=0; posi<confIndex.numUndefined; posi++) {
            int pos = confIndex.undefinedPos[posi];
            undefinedOrder.add(pos);
            scores.put(pos, scorePos(confIndex, rcs, pos));
        }

        List<Integer> breakList = calcStrandBreaks(undefinedOrder, confspace);

        for ( int breakIndex=0; breakIndex < breakList.size()-1 ; breakIndex++ ) {
            // sort positions in order of decreasing score
            Collections.sort(undefinedOrder.subList(breakList.get(breakIndex), breakList.get(breakIndex+1)), new Comparator<Integer>() {

                @Override
                public int compare(Integer pos1, Integer pos2) {
                    double score1 = scores.get(pos1);
                    double score2 = scores.get(pos2);
                    // NOTE: use reverse order for decreasing sort
                    return Double.compare(score2, score1);
                }
            });
        }

        // prepend the defined positions to build the full order
        List<Integer> order = new ArrayList<>();
        for (int posi=0; posi<confIndex.numDefined; posi++) {
            int pos = confIndex.definedPos[posi];
            order.add(pos);
        }
        order.addAll(undefinedOrder);


        return order;
    }

    private List<Integer> calcStrandBreaks(List<Integer> indexList, SimpleConfSpace confspace){
        /** Computes the indices at which positions change strands
         *
         */
        List strandList = indexList.stream()
                .map((indx) -> confspace.positions.get(indx).strand)
                .collect(Collectors.toList());
        List<Integer> breakList = new ArrayList<Integer>();
        breakList.add(0);
        for (int i = 0; i < strandList.size()-1; i++) {
            if ( strandList.get(i) != strandList.get(i+1) ){
               breakList.add(i+1);
            }
        }
        // add integers when strand switches
        breakList.add(indexList.size());

        return breakList;
    }
}
