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
	public void compute(long maxNumConfs) {
		if (!status.canContinue()) {
			throw new IllegalStateException("can't continue from status " + status);
		}

		// avoid overflow by capping min score
		ScoredConf conf; final double minScore = -4000;
		if ((conf = scoreConfs.next()) != null) {
			double score = conf.getScore();
			if(score < minScore) score = minScore;
			values.qstar = boltzmann.calc(score);
			//values.qstar = BigDecimalUtil.exp(BigDecimal.valueOf(-conf.getScore()/PoissonBoltzmannEnergy.constRT), BigDecimalUtil.getScale());
		}

		this.targetEpsilon = 0;
		status = Status.Estimated;
	}

}
