package edu.duke.cs.osprey.astar.conf.order;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

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
		return posOrder.get(confIndex.getNode().getLevel());
	}
	
	private List<Integer> calcPosOrder(ConfIndex confIndex, RCs rcs) {
		
		// init permutation array
		List<Integer> order = new ArrayList<>();
		for (int pos=0; pos<rcs.getNumPos(); pos++) {
			order.add(pos);
		}
		
		// calculate all the scores
		List<Double> scores = new ArrayList<>();
		for (int pos=0; pos<rcs.getNumPos(); pos++) {
			scores.add(scorePos(confIndex, rcs, pos));
		}
		
		// sort positions in order of decreasing score
		Collections.sort(order, new Comparator<Integer>() {

			@Override
			public int compare(Integer pos1, Integer pos2) {
				double score1 = scores.get(pos1);
				double score2 = scores.get(pos2);
				// NOTE: use reverse order for decreasing sort
				return Double.compare(score2, score1);
			}
		});
			
		return order;
	}

	double scorePos(ConfIndex confIndex, RCs rcs, int pos) {
		
		// check all the RCs at this pos and aggregate the energies
		double parentScore = confIndex.getNode().getScore();
		double reciprocalSum = 0;
		for (int rc : rcs.get(pos)) {
			double childScore = gscorer.calcDifferential(confIndex, rcs, pos, rc)
				+ hscorer.calcDifferential(confIndex, rcs, pos, rc);
			reciprocalSum += 1.0/(childScore - parentScore);
		}
		return 1.0/reciprocalSum;
	}
}
