package edu.duke.cs.osprey.astar.conf.order;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;

public class DynamicHMeanAStarOrder implements AStarOrder {
	
	private AStarScorer gscorer;
	private AStarScorer hscorer;
	
	@Override
	public void setScorers(AStarScorer gscorer, AStarScorer hscorer) {
		this.gscorer = gscorer;
		this.hscorer = hscorer;
	}

	@Override
	public int getNextPos(ConfIndex confIndex, RCs rcs) {
		
		int bestPos = -1;
		double bestScore = Double.NEGATIVE_INFINITY;
		
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			
			int pos = confIndex.getUndefinedPos()[i];
			double score = scorePos(confIndex, rcs, pos);
			
			if (score > bestScore) {
				bestScore = score;
				bestPos = pos;
			}
		}
			
		if (bestPos == -1) {
			throw new RuntimeException("ERROR: Can't find next position for dynamic A*");
		}
		
		return bestPos;
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
