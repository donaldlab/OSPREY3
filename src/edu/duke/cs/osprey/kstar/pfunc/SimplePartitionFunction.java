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

package edu.duke.cs.osprey.kstar.pfunc;

import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.JvmMem;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.math.BigDecimal;
import java.math.BigInteger;


/**
 * This partition function calculator estimates the partition function value
 * by alternating between two operations:
 * 	1. compute lower energy bounds on conformations to refine the upper pfunc bound
 * and:
 * 	2. compute upper energy bounds on conformations to refine the lower pfunc bound
 *
 * This pfunc implementation prefers to do operation 1 only when the thread pool is busy
 * with operation 2. This approach works really well when operation 1 is orders of
 * magnitude less expensive than operation 2 (the typical case), but doesn't perform
 * well when operation 2 becomes much faster (say if we e.g. were reading energies out
 * of a cache). An attempt to solve this problem is implemented in the
 * SimplerPartitionFunction class.
 *
 * Operation 2 is performed while the thread pool is busy, but will eventually
 * cease even while the thread pool is busy if we determine that the upper bound is
 * already tight enough.
 */
public class SimplePartitionFunction implements PartitionFunction {

    public final ConfEnergyCalculator ecalc;

    /**
	 * number of conformations to score per batch when refining the partition function upper bound
	 *
	 * for single-threaded parallelism, this setting determines how many conformations we should score
	 * to refine the upper bound *per each* conformation that is minimized to refine the lower bound.
	 *
	 * for multi-threaded (including GPUs) parallelism, this setting is less important, since the
	 * upper bound calculator will automatically run on the main thread during the time the main
	 * thread would have normally waited for the energy minimizer thread(s).
	 **/
    public int scoreConfsBatchSize = 1000;

	private double targetEpsilon = Double.NaN;
	private BigDecimal stabilityThreshold = BigDecimal.ZERO;
	private Status status = null;
	private Values values = null;
	private LowerBoundCalculator lowerBound;
	private UpperBoundCalculator upperBound;
	private ConfListener confListener = null;
	private boolean isReportingProgress = false;
	private Stopwatch stopwatch = new Stopwatch().start();
	private ConfDB confDB = null;


	public SimplePartitionFunction(ConfEnergyCalculator ecalc) {
		this.ecalc = ecalc;
	}
	
	@Override
	public void setReportProgress(boolean val) {
		isReportingProgress = val;
	}
	
	@Override
	public void setConfListener(ConfListener val) {
		confListener = val;
	}
	
	@Override
	public Status getStatus() {
		return status;
	}
	
	@Override
	public Values getValues() {
		return values;
	}
	
	@Override
	public int getNumConfsEvaluated() {
		return lowerBound.numConfsEnergied;
	}
	
	@Override
	public int getParallelism() {
		return ecalc.tasks.getParallelism();
	}

	@Override
	public void init(ConfSearch confSearch, BigInteger numConfsBeforePruning, double targetEpsilon) {

		this.targetEpsilon = targetEpsilon;

		status = Status.Estimating;
		values = Values.makeFullRange();

		// don't explicitly check the pruned confs, just lump them together with the un-enumerated confs
		values.pstar = BigDecimal.ZERO;

		// split the confs between the bound calculators
		ConfSearch.MultiSplitter confsSplitter = new ConfSearch.MultiSplitter(confSearch);
		lowerBound = new LowerBoundCalculator(confsSplitter.makeStream(), ecalc);
		lowerBound.confDB = confDB;
		upperBound = new UpperBoundCalculator(confsSplitter.makeStream(), numConfsBeforePruning);
	}

	@Override
	public void setStabilityThreshold(BigDecimal val) {
		this.stabilityThreshold = val;
	}

	@Override
	public void compute(int maxNumConfs) {

		if (status == null) {
			throw new IllegalStateException("pfunc was not initialized. Call init() before compute()");
		}

		// WARNING: if this is too big, we'll never reach the pfunc epsilon!
		final double upperBoundEpsilon = 0.00001;

		// step the upper bound calculator at least once,
		// so we try to get a non-inf upper bound before working on the lower bound
		upperBound.scoreNextConf();
		values.qprime = upperBound.totalBound.subtract(lowerBound.weightedScoreSum);

		int initialNumConfsScored = lowerBound.numConfsScored;

		while (true) {

			// don't race the calculators while we're checking stopping criteria
			synchronized (this) {

				// did we drop below the stability threshold?
				if (lowerBound.numConfsEnergied > 0 && stabilityThreshold != null && MathTools.isLessThan(values.calcUpperBound(), stabilityThreshold)) {
					status = Status.Unstable;
					break;
				}

				// did we hit the epsilon target?
				if (values.getEffectiveEpsilon() <= targetEpsilon) {
					status = Status.Estimated;
					break;
				}

				// should we stop anyway?
				if (!status.canContinue() || lowerBound.numConfsScored - initialNumConfsScored >= maxNumConfs) {
					break;
				}
			}

			// nope, need to do some more work

			if (ecalc.tasks instanceof ThreadPoolTaskExecutor) {

				// while we're waiting on energy threads, refine q' a bit on the main thread
				while (upperBound.delta > upperBoundEpsilon && ecalc.tasks.isBusy()) {
					upperBound.run(scoreConfsBatchSize, upperBoundEpsilon);
				}

			} else {

				// we're stuck with just the main thread
				// meaning we need to refine q' explicitly before computing energies
				upperBound.run(scoreConfsBatchSize, upperBoundEpsilon);
			}

			// refine the lower bound
			LowerBoundCalculator.Status lbStatus = lowerBound.energyNextConfAsync((econf) -> {

				// got an energy for the conf

				// we might be on the listener thread, so lock to keep from racing the main thread
				synchronized (SimplePartitionFunction.this) {

					// update pfunc values
					values.qstar = lowerBound.weightedEnergySum;
					values.qprime = upperBound.totalBound.subtract(lowerBound.weightedScoreSum);
				}

				// report progress if needed
				if (isReportingProgress) {
					System.out.println(String.format("conf:%4d, score:%12.6f, energy:%12.6f, q*:%12e, q':%12e, epsilon:%.6f, time:%10s, heapMem:%s, extMem:%s",
						lowerBound.numConfsEnergied, econf.getScore(), econf.getEnergy(), values.qstar, values.qprime, values.getEffectiveEpsilon(),
						stopwatch.getTime(2),
						JvmMem.getOldPool(),
						ExternalMemory.getUsageReport()
					));
				}

				// report confs if needed
				if (confListener != null) {
					confListener.onConf(econf);
				}
			});

			// see if the lower bound calculator stalled and wait if it did
			switch (lbStatus) {
				case OutOfConfs:
					status = Status.OutOfConformations;
					lowerBound.waitForFinish();
				break;
				case OutOfLowEnergies:
					status = Status.OutOfLowEnergies;
					lowerBound.waitForFinish();
				break;
			}

			// make sure the lower bound calc doesn't get ahead of the upper bound calc
			// otherwise, the bounds will be wrong
			while (lowerBound.numConfsScored > upperBound.numScoredConfs) {
				int oldNum = upperBound.numScoredConfs;
				upperBound.scoreNextConf();

				// make sure it's even possible for the upper bound calculator to score more confs
				if (upperBound.numScoredConfs <= oldNum) {
					throw new Error(String.format(
						"The lower bound calculator apparently scored more conformations (%d) than is possible"
						+ " for the upper bound calculator (%s). This is definitely a bug",
						lowerBound.numConfsScored,
						upperBound.tree.getNumConformations().toString()
					));
				}
			}
		}

		// if we're still waiting on energy threads, refine q' a bit more on the main thread
		if (ecalc.tasks instanceof ThreadPoolTaskExecutor) {
			while (upperBound.delta > upperBoundEpsilon && ecalc.tasks.isWorking()) {
				upperBound.run(scoreConfsBatchSize, upperBoundEpsilon);
			}
		}

		// wait for any remaining async minimizations to finish
		ecalc.tasks.waitForFinish();

		// after we've decided this pfunc is sufficiently estimated,
		// sometimes extra energies can arrive asynchronously
		// and push us into unstable territory,
		// so check for that after waiting on the ecalc to finish
		if (status != Status.Unstable && stabilityThreshold != null && MathTools.isLessThan(values.calcUpperBound(), stabilityThreshold)) {
			status = Status.Unstable;
		}
	}
}
