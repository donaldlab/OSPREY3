package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.TuplesIndex;


public class LUTEGScorer implements AStarScorer {

	public final LUTEConfEnergyCalculator ecalc;

	public LUTEGScorer(LUTEConfEnergyCalculator ecalc) {
		this.ecalc = ecalc;
	}

	@Override
	public AStarScorer make() {
		return new LUTEGScorer(ecalc);
	}

	@Override
	public double calc(ConfIndex confIndex, RCs rcs) {

		/* convert the conf index into a conf
			yeah, we'll take a slight performance hit doing this,
			but the higher-order tuple indices are optimized around fast pos->rc lookups,
			but ConfIndex optimizes fast singles/pairs enumeration instead
		*/
		int[] conf = Conf.make(confIndex);

		try {
			return ecalc.calcEnergy(conf);
		} catch (TuplesIndex.NoSuchTupleException ex) {
			// conf has a pruned tuple, can't score it
			return Double.POSITIVE_INFINITY;
		}
	}
}
