package edu.duke.cs.osprey.astar.conf.order;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;

public class StaticScoreHMeanAStarOrder implements AStarOrder {
	
	private AStarScorer gscorer;
	private AStarScorer hscorer;
	private List<Integer> posOrder;
	
	public StaticScoreHMeanAStarOrder() {
		posOrder = null;
	}
	
	@Override
	public void setScorers(AStarScorer gscorer, AStarScorer hscorer) {
		this.gscorer = gscorer;
		this.hscorer = hscorer;
	}

	@Override
	public int getNextPos(ConfIndex confIndex, RCs rcs) {
		if (posOrder == null) {
			posOrder = calcPosOrder(confIndex, rcs);
		}
		return posOrder.get(confIndex.node.getLevel());
	}
	
	private List<Integer> calcPosOrder(ConfIndex confIndex, RCs rcs) {
		
		// init permutation array with only undefined positions and score them
		List<Integer> undefinedOrder = new ArrayList<Integer>();
		Map<Integer,Double> scores = new TreeMap<>();
		for (int posi=0; posi<confIndex.numUndefined; posi++) {
			int pos = confIndex.undefinedPos[posi];
			undefinedOrder.add(pos);
			scores.put(pos, scorePos(confIndex, rcs, pos));
		}
		
		// sort positions in order of decreasing score
		Collections.sort(undefinedOrder, new Comparator<Integer>() {

			@Override
			public int compare(Integer pos1, Integer pos2) {
				double score1 = scores.get(pos1);
				double score2 = scores.get(pos2);
				// NOTE: use reverse order for decreasing sort
				return Double.compare(score2, score1);
			}
		});
		
		// prepend the defined positions to build the full order
		List<Integer> order = new ArrayList<>();
		for (int posi=0; posi<confIndex.numDefined; posi++) {
			int pos = confIndex.definedPos[posi];
			order.add(pos);
		}
		order.addAll(undefinedOrder);
		
		return order;
	}

	double scorePos(ConfIndex confIndex, RCs rcs, int pos) {
		
		// check all the RCs at this pos and aggregate the energies
		double parentScore = confIndex.node.getScore();
		double reciprocalSum = 0;
		for (int rc : rcs.get(pos)) {
			double childScore = gscorer.calcDifferential(confIndex, rcs, pos, rc)
				+ hscorer.calcDifferential(confIndex, rcs, pos, rc);
			reciprocalSum += 1.0/(childScore - parentScore);
		}
		return 1.0/reciprocalSum;
	}
}
