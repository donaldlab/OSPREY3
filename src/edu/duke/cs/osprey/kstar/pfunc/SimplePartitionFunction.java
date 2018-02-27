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

	public final ConfSearch confSearch;
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
	private BigDecimal upperBoundThreshold = null;
	private Status status = null;
	private Values values = null;
	private LowerBoundCalculator lowerBound;
	private UpperBoundCalculator upperBound;
	private ConfListener confListener = null;
	private boolean isReportingProgress = false;
	private Stopwatch stopwatch = new Stopwatch().start();
	private ConfDB confDB = null;


	public SimplePartitionFunction(ConfSearch confSearch, ConfEnergyCalculator ecalc) {
		this.confSearch = confSearch;
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
	public void init(double targetEpsilon) {
		init(targetEpsilon, BigDecimal.ZERO);
	}

	@Override
	public void init(double targetEpsilon, BigDecimal stabilityThreshold) {

		this.targetEpsilon = targetEpsilon;
		this.upperBoundThreshold = stabilityThreshold;
		
		status = Status.Estimating;
		values = Values.makeFullRange();

		// NOTE: don't use DEE with this pfunc calculator
		// DEE actually makes the problem harder to solve, not easier
		// because then we have to deal with p*
		// if we just don't prune with DEE, the usual q' calculator will handle those confs that would have been pruned
		// not using DEE won't really slow us down either, since our A* is fast enough without it
		values.pstar = BigDecimal.ZERO;

		// split the confs between the bound calculators
		ConfSearch.Splitter confsSplitter = new ConfSearch.Splitter(confSearch);
		lowerBound = new LowerBoundCalculator(confsSplitter.makeStream(), ecalc);
		lowerBound.confDB = confDB;
		upperBound = new UpperBoundCalculator(confsSplitter.makeStream());
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
				if (lowerBound.numConfsEnergied > 0 && upperBoundThreshold != null && MathTools.isLessThan(values.calcUpperBound(), upperBoundThreshold)) {
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
						confSearch.getNumConformations().toString()
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
		if (status != Status.Unstable && upperBoundThreshold != null && MathTools.isLessThan(values.calcUpperBound(), upperBoundThreshold)) {
			status = Status.Unstable;
		}
	}
}
