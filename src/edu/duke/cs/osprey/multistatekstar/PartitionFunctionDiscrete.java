package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.control.ConfEnergyCalculator.Async;

import java.math.BigDecimal;
import java.math.BigInteger;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class PartitionFunctionDiscrete extends PartitionFunctionMinimized {

	public static boolean DEBUG = false;

	public PartitionFunctionDiscrete(
			EnergyMatrix emat, 
			PruningMatrix pmat, 
			PruningMatrix invmat, 
			ConfSearchFactory confSearchFactory,
			Async ecalc
			) {
		super(emat, pmat, invmat, confSearchFactory, ecalc);
	}
	
	private ScoredConf findRigidGMEC() {
		return energyConfs.next();
	}
	
	@Override
	public void init(double targetEpsilon) {
		super.init(targetEpsilon);
		scoreConfs = null;
		
		//only applies to the final object
		if(computeGMECRatio) {
			ScoredConf conf = findRigidGMEC();
			double gmecEnergy = conf == null ? Double.POSITIVE_INFINITY : conf.getScore();
			values.qstar = boltzmann.calc(gmecEnergy);
			numConfsEvaluated++;
			if (confListener != null) {
				confListener.onConf(conf);
			}
			status = Status.Estimated;
		}
	}

	@Override
	public void compute(int maxNumConfs) {

		if (!status.canContinue()) {
			throw new IllegalStateException("can't continue from status " + status);
		}

		ScoredConf conf;
		BigDecimal scoreWeight;
		
		int stopAtConf = numConfsEvaluated + maxNumConfs;

		double lastScore = Double.NEGATIVE_INFINITY;

		while (true) {

			// should we keep going?
			if (!status.canContinue() || numConfsEvaluated >= stopAtConf) {
				break;
			}

			if ((conf = energyConfs.next()) == null) {
				if(status != Status.Estimated) status = Status.NotEnoughConformations;
				break;
			}

			numConfsEvaluated++;

			if(DEBUG) {
				if(conf.getScore() < lastScore)
					throw new RuntimeException("ERROR: scores must be non-decreasing: score: "+conf.getScore()+ ", lastScore: "+lastScore);
				lastScore = conf.getScore();
			}

			scoreWeight = boltzmann.calc(conf.getScore());
			numConfsToScore = numConfsToScore.subtract(BigInteger.ONE);

			if (scoreWeight.compareTo(BigDecimal.ZERO) == 0) {
				values.qprime = updateQprime(scoreWeight);
				double effectiveEpsilon = getEffectiveEpsilon();	
				if (!Double.isNaN(effectiveEpsilon) && effectiveEpsilon <= targetEpsilon) status = Status.Estimated;
				else if(status != Status.Estimated) status = Status.NotEnoughFiniteEnergies;
				break;
			}

			if(status == Status.Estimating) {

				values.qstar = values.qstar.add(scoreWeight);
				values.qprime = updateQprime(scoreWeight);

				// report progress if needed
				if (isReportingProgress && numConfsEvaluated % 1024 == 0) {
					phase1Output(conf);
				}

				// report confs if needed
				if (confListener != null) {
					confListener.onConf(conf);
				}

				// update status if needed
				double effectiveEpsilon = getEffectiveEpsilon();
				if(Double.isNaN(effectiveEpsilon)) {
					status = Status.NotEnoughFiniteEnergies;
				}
				else if (effectiveEpsilon <= targetEpsilon) {
					status = Status.Estimated;
					if (isReportingProgress) phase1Output(conf);//just to let the user know we reached epsilon
				}
			}
		}
	}

	@Override
	public void compute(BigDecimal targetScoreWeights, int maxNumConfs) {

		if (!status.canContinue()) {
			throw new IllegalStateException("can't continue from status " + status);
		}

		ScoredConf conf;
		BigDecimal scoreWeight;

		int stopAtConf = numConfsEvaluated + maxNumConfs;
		while (true) {

			// should we keep going?
			if (!status.canContinue() 
					|| qstarScoreWeights.compareTo(targetScoreWeights) >= 0
					|| numConfsEvaluated >= stopAtConf) {
				break;
			}

			if ((conf = energyConfs.next()) == null) {
				if(status != Status.Estimated) status = Status.NotEnoughConformations;
				break;
			}

			numConfsEvaluated++;

			scoreWeight = boltzmann.calc(conf.getScore());
			numConfsToScore = numConfsToScore.subtract(BigInteger.ONE);

			if (scoreWeight.compareTo(BigDecimal.ZERO) == 0) {
				values.qprime = updateQprime(scoreWeight);
				double effectiveEpsilon = getEffectiveEpsilon();	
				if (!Double.isNaN(effectiveEpsilon) && effectiveEpsilon <= targetEpsilon) status = Status.Estimated;
				else if(status != Status.Estimated) status = Status.NotEnoughFiniteEnergies;
				break;
			}

			if(status == Status.Estimating) {
				// get the boltzmann weight
				qstarScoreWeights = qstarScoreWeights.add(scoreWeight);	

				// update pfunc state
				values.qstar = values.qstar.add(scoreWeight);
				values.qprime = updateQprime(scoreWeight);
				BigDecimal pdiff = targetScoreWeights.subtract(qstarScoreWeights);

				// report progress if needed
				if (isReportingProgress && numConfsEvaluated % 1024 == 0) {
					phase2Output(conf, pdiff);
				}

				// report confs if needed
				if (confListener != null) {
					confListener.onConf(conf);
				}

				// update status if needed
				double effectiveEpsilon = getEffectiveEpsilon();
				if(Double.isNaN(effectiveEpsilon)) {
					status = Status.NotEnoughFiniteEnergies;
				}
				else if (effectiveEpsilon <= targetEpsilon) {
					status = Status.Estimated;
					if (isReportingProgress) phase2Output(conf, pdiff);
				}
			}
		}
	}

	protected BigDecimal updateQprime(BigDecimal val) {
		return val.multiply(new BigDecimal(numConfsToScore.toString()));
	}

	protected double getEffectiveEpsilon() {
		return values.getEffectiveEpsilon();
	}
}
