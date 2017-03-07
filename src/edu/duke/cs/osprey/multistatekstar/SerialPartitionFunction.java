package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.control.ConfEnergyCalculator.Async;

import java.math.BigDecimal;
import java.math.BigInteger;

import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class SerialPartitionFunction extends ParallelPartitionFunction {

	public SerialPartitionFunction(EnergyMatrix emat, PruningMatrix pmat, ConfSearchFactory confSearchFactory,
			Async ecalc) {
		super(emat, pmat, confSearchFactory, ecalc);
	}

	@Override
	public void init(double targetEpsilon) {

	}

	@Override
	public void compute(int maxNumConfs) {

		if (!status.canContinue()) {
			throw new IllegalStateException("can't continue from status " + status);
		}

		int stopAtConf = numConfsEvaluated + maxNumConfs;
		while (true) {

			ScoredConf conf;
			BigDecimal scoreWeight;

			// should we keep going?
			if (!status.canContinue() || numConfsEvaluated >= stopAtConf) {
				break;
			}

			if ((conf = energyConfs.next()) == null) {
				if(status != Status.Estimated) status = Status.NotEnoughConformations;
				break;
			}

			numConfsEvaluated++;

			scoreWeight = boltzmann.calc(conf.getScore());

			if (scoreWeight.compareTo(BigDecimal.ZERO) == 0) {
				if(status != Status.Estimated) status = Status.NotEnoughFiniteEnergies;
				break;
			}

			if(status == Status.Estimating) {

				numConfsToScore = numConfsToScore.subtract(BigInteger.ONE);
				
				values.qstar = values.qstar.add(scoreWeight);
				values.qprime = updateQprime(scoreWeight);

				// report progress if needed
				if (isReportingProgress && numConfsEvaluated % ecalc.getParallelism() == 0) {
					phase1Output(conf);
				}
				
				// report confs if needed
				if (confListener != null) {
					confListener.onConf(conf);
				}

				// update status if needed
				if (values.getEffectiveEpsilon() <= targetEpsilon) {
					status = Status.Estimated;
					phase1Output(conf);//just to let the user know we reached epsilon
				}
			}
		}
	}

	protected BigDecimal updateQprime(BigDecimal val) {
		return val.multiply(new BigDecimal(numConfsToScore.toString()));
	}
}
