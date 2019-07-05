/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.math.BigInteger;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.gmec.ConfSearchFactory;
import edu.duke.cs.osprey.gmec.GMECConfEnergyCalculator.Async;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class PartitionFunctionDiscrete extends PartitionFunctionMinimized {

	public PartitionFunctionDiscrete(
			EnergyMatrix emat, 
			PruningMatrix pmat, 
			PruningMatrix invmat, 
			ConfSearchFactory confSearchFactory,
			Async ecalc
			) {
		super(emat, pmat, invmat, confSearchFactory, ecalc);
	}

	@Override
	public void init(ConfSearch confSearch, BigInteger numConfsBeforePruning, double targetEpsilon) {
		super.init(confSearch, numConfsBeforePruning, targetEpsilon);
		energyConfs = null;
	}

	@Override
	public void compute(int maxNumConfs) {

		if (!status.canContinue()) {
			throw new IllegalStateException("can't continue from status " + status);
		}

		ScoredConf conf;
		BigDecimal scoreWeight;
		int stopAtConf = numConfsEvaluated + maxNumConfs;

		while (true) {

			// should we keep going?
			if (!status.canContinue() || numConfsEvaluated >= stopAtConf) {
				break;
			}

			if ((conf = scoreConfs.nextConf()) == null) {
				if(status != Status.Estimated) status = Status.OutOfConformations;
				break;
			}

			numConfsEvaluated++;

			scoreWeight = boltzmann.calc(conf.getScore());

			if (scoreWeight.compareTo(BigDecimal.ZERO) == 0) {
				if(status != Status.Estimated) status = Status.OutOfLowEnergies;
				break;
			}

			if(status == Status.Estimating) {

				numConfsToScore = numConfsToScore.subtract(BigInteger.ONE);

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
				if (getEffectiveEpsilon() <= targetEpsilon) {
					status = Status.Estimated;
					phase1Output(conf);//just to let the user know we reached epsilon
				}
			}
		}
	}

	public void compute(BigDecimal targetScoreWeights) {

		if (!status.canContinue()) {
			throw new IllegalStateException("can't continue from status " + status);
		}

		ScoredConf conf;
		BigDecimal scoreWeight;

		while (true) {

			// should we keep going?
			if (!status.canContinue() || qstarScoreWeights.compareTo(targetScoreWeights) >= 0) {
				break;
			}

			if ((conf = scoreConfs.nextConf()) == null) {
				if(status != Status.Estimated) status = Status.OutOfConformations;
				break;
			}

			numConfsEvaluated++;

			scoreWeight = boltzmann.calc(conf.getScore());

			if (scoreWeight.compareTo(BigDecimal.ZERO) == 0) {
				if(status != Status.Estimated) status = Status.OutOfLowEnergies;
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
				if (getEffectiveEpsilon() <= targetEpsilon) {
					status = Status.Estimated;
					phase2Output(conf, pdiff);
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
