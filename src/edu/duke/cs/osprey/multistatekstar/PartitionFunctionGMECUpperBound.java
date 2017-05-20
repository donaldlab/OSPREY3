package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.control.ConfEnergyCalculator.Async;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class PartitionFunctionGMECUpperBound extends PartitionFunctionDiscrete {

	public PartitionFunctionGMECUpperBound(
			EnergyMatrix emat, 
			PruningMatrix pmat,
			PruningMatrix invmat,
			ConfSearchFactory confSearchFactory, 
			Async ecalc
			) {
		super(emat, pmat, invmat, confSearchFactory, ecalc);
	}

	@Override
	public void init(double targetEpsilon) {
		this.targetEpsilon = targetEpsilon;

		status = Status.Estimating;
		values = new Values();
		
		// make the search tree for computing q*
		ConfSearch tree = confSearchFactory.make(emat, pmat);
		if(tree instanceof ConfAStarTree) ((ConfAStarTree)tree).stopProgress();
		ConfSearch.Splitter confsSplitter = new ConfSearch.Splitter(tree);
		scoreConfs = confsSplitter.makeStream();
	}

	@Override
	public void compute(int maxNumConfs) {
		if (!status.canContinue()) {
			throw new IllegalStateException("can't continue from status " + status);
		}

		ScoredConf conf;
		if ((conf = scoreConfs.next()) != null) {
			values.qstar = boltzmann.calc(conf.getScore());
		}

		this.targetEpsilon = 0;
		status = Status.Estimated;
	}

}
