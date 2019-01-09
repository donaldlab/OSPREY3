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

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.tools.*;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Supplier;


/**
 * This partition function calculator estimates the partition function value
 * by alternating between two operations:
 * 	1. compute lower energy bounds on conformations to refine the pfunc upper bound
 * and:
 * 	2. compute upper energy bounds on conformations to refine the pfunc lower bound
 *
 * This implementation will always do operation 1 first, then operation 2 to get off
 * the initial "flat spot" the pfunc surface.
 *
 * After that, this implementation is an attempt at a gradient-descent type partition
 * function calculator. At each step, the operation that maximizes the drop in delta
 * (according to slope sampling heuristics) is chosen.
 *
 * This implementation should also perform much better whe operation 2 is
 * not orders of magnitude slower than operation 1 (when e.g. we're reading
 * energies out of a cache).
 */
public class GradientDescentPfunc implements PartitionFunction.WithConfTable, PartitionFunction.WithExternalMemory {

	private static class State {

		BigDecimal numConfs;

		// upper bound (score axis) vars
		long numScoredConfs = 0;
		BigDecimal upperScoreWeightSum = BigDecimal.ZERO;
		BigDecimal minUpperScoreWeight = MathTools.BigPositiveInfinity;
		double minUpperScore = Double.NEGATIVE_INFINITY;

		// lower bound (energy axis) vars
		long numEnergiedConfs = 0;
		BigDecimal lowerScoreWeightSum = BigDecimal.ZERO;
		BigDecimal energyWeightSum = BigDecimal.ZERO;
		BigDecimal minLowerScoreWeight = MathTools.BigPositiveInfinity;
		BigDecimal cumulativeZReduction = BigDecimal.ZERO;
		ArrayList<Integer> minList = new ArrayList<Integer>();
		BigDecimal firstScoreWeight = BigDecimal.ZERO;

		// estimate of inital rates
		// (values here aren't super imporant since they get tuned during execution,
		// but make scores much faster than energies so we don't get a slow start)
		double scoreOps = 100.0;
		double energyOps = 1.0;

		double prevDelta = 1.0;
		double dEnergy = -1.0;
		double dScore = -1.0;

		State(BigInteger numConfs) {
			this.numConfs = new BigDecimal(numConfs);
		}

		double calcDelta() {
			BigDecimal upperBound = getUpperBound();
			if (MathTools.isZero(upperBound) || MathTools.isInf(upperBound)) {
				return 1.0;
			}
			return new BigMath(PartitionFunction.decimalPrecision)
				.set(upperBound)
				.sub(getLowerBound())
				.div(upperBound)
				.get()
				.doubleValue();
		}

		public BigDecimal getLowerBound() {
			return energyWeightSum;
		}

		public void printBoundStats() {
            System.out.println("Num confs: " + String.format("%12e",numConfs));
            System.out.println("Num Scored confs: " + String.format("%4d",numScoredConfs));
            String upperScoreString = minUpperScoreWeight.toString();
            String upperSumString = upperScoreWeightSum.toString();
            if(!MathTools.isInf(minUpperScoreWeight))
                upperScoreString = String.format("%12e",minUpperScoreWeight);
            if(!MathTools.isInf(upperScoreWeightSum))
                upperSumString = String.format("%12e",upperScoreWeightSum);
            System.out.println("Conf bound: " + upperScoreString);
            System.out.println("Scored weight bound:"+ upperSumString);
		}

		public BigDecimal getUpperBound() {

			return new BigMath(PartitionFunction.decimalPrecision)

				// unscored bound
				.set(numConfs)
				.sub(numScoredConfs)
				.mult(minUpperScoreWeight)

				// with scored bound
				.add(upperScoreWeightSum)

				// but replace weights that have energies
				.sub(lowerScoreWeightSum)
				.add(energyWeightSum)

				.get();
		}
		public BigDecimal getUpperBoundNoE() {

			return new BigMath(PartitionFunction.decimalPrecision)

					// unscored bound
					.set(numConfs)
					.sub(numScoredConfs)
					.mult(minUpperScoreWeight)

					// with scored bound
					.add(upperScoreWeightSum)

					.get();
		}

		boolean epsilonReached(double targetEpsilon) {
			return calcDelta() <= targetEpsilon;
		}

		boolean isStable(BigDecimal stabilityThreshold) {
			return numEnergiedConfs <= 0 || stabilityThreshold == null || MathTools.isGreaterThanOrEqual(getUpperBound(), stabilityThreshold);
		}

		boolean hasLowEnergies() {
			return MathTools.isGreaterThan(minLowerScoreWeight,  BigDecimal.ZERO);
		}

		@Override
		public String toString() {
			return String.format("upper: count %d  sum %e  min %e     lower: count %d  score sum %e  energy sum %e",
				numScoredConfs, upperScoreWeightSum, minUpperScoreWeight,
				numEnergiedConfs, lowerScoreWeightSum, energyWeightSum
			);
		}
	}

	private static enum Step {
		None,
		Score,
		Energy
	}


    public final ConfEnergyCalculator ecalc;

	private double targetEpsilon = Double.NaN;
	private BigDecimal stabilityThreshold = BigDecimal.ZERO;
	private ConfListener confListener = null;
	private boolean isReportingProgress = false;
	private Stopwatch stopwatch = new Stopwatch().start();
	private ConfSearch scoreConfs = null;
	private ConfSearch energyConfs = null;
	private BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

	private Status status = null;
	private Values values = null;
	private State state = null;

	private boolean hasEnergyConfs = true;
	private boolean hasScoreConfs = true;
	private long numEnergyConfsEnumerated = 0;
	private long numScoreConfsEnumerated = 0;

	private ConfDB.ConfTable confTable = null;

	private boolean useExternalMemory = false;
	private RCs rcs = null;

	private PfuncSurface surf = null;
	private PfuncSurface.Trace trace = null;

	public GradientDescentPfunc(ConfEnergyCalculator ecalc) {
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
		// TODO: this might overflow for big pfunc calculations, upgrade the interface type?
		return (int)state.numEnergiedConfs;
	}

	public int getNumConfsScored() {
		return (int) state.numScoredConfs;
	}

	@Override
	public int getParallelism() {
		return ecalc.tasks.getParallelism();
	}

	@Override
	public void setConfTable(ConfDB.ConfTable val) {
		confTable = val;
	}

	@Override
	public void setUseExternalMemory(boolean val, RCs rcs) {
		this.useExternalMemory = val;
		this.rcs = rcs;
	}

	public void traceTo(PfuncSurface val) {
		surf = val;
	}

	@Override
	public void init(ConfSearch confSearch, BigInteger numConfsBeforePruning, double targetEpsilon) {

		init(numConfsBeforePruning, targetEpsilon);

		// split the confs between the upper and lower bounds
		ConfSearch.Splitter confsSplitter = new ConfSearch.Splitter(confSearch, useExternalMemory, rcs);
		scoreConfs = confsSplitter.first;
		energyConfs = confsSplitter.second;
	}

	@Override
	public void init(ConfSearch upperBoundConfs, ConfSearch lowerBoundConfs, BigInteger numConfsBeforePruning, double targetEpsilon) {

		init(numConfsBeforePruning, targetEpsilon);

		this.scoreConfs = upperBoundConfs;
		this.energyConfs = lowerBoundConfs;
	}

	private void init(BigInteger numConfsBeforePruning, double targetEpsilon) {

		if (targetEpsilon <= 0.0) {
			throw new IllegalArgumentException("target epsilon must be greater than zero");
		}

		this.targetEpsilon = targetEpsilon;

		// init state
		status = Status.Estimating;
		state = new State(numConfsBeforePruning);
		values = Values.makeFullRange();
		// don't explicitly check the pruned confs, just lump them together with the un-enumerated confs
		values.pstar = BigDecimal.ZERO;

		hasEnergyConfs = true;
		hasScoreConfs = true;
		numEnergyConfsEnumerated = 0;
		numScoreConfsEnumerated = 0;
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
		if (!status.canContinue()) {
			return;
		}

		// start a trace if needed
		if (surf != null) {
			trace = surf.new Trace();
		}

		boolean keepStepping = true;
		for (int numConfsEnergied=0; numConfsEnergied<maxNumConfs; /* don't increment here */) {

			// which way should we step, and how far?
			Step step = Step.None;
			int numScores = 0;
			synchronized (this) { // don't race the listener thread

				// should we even keep stepping?
				keepStepping = keepStepping
					&& !state.epsilonReached(targetEpsilon)
					&& state.isStable(stabilityThreshold)
					&& state.hasLowEnergies();
				if (!keepStepping) {
					break;
				}

				// just in case...
				if (Double.isNaN(state.dEnergy) || Double.isNaN(state.dScore)) {
					throw new Error("Can't determine gradient of delta surface. This is a bug.");
				}

				boolean scoreAheadOfEnergy = numEnergyConfsEnumerated < numScoreConfsEnumerated;
				boolean energySteeperThanScore = state.dEnergy <= state.dScore;

				if (hasEnergyConfs && ((scoreAheadOfEnergy && energySteeperThanScore) || !hasScoreConfs)) {

					step = Step.Energy;

				} else if (hasScoreConfs) {

					step = Step.Score;

					// how many scores should we weight?
					// (target a similar amount of time as energy calculation, but at least 10 ms)
					double scoringSeconds = Math.max(0.1/state.energyOps, 0.01);
					numScores = Math.max((int)(scoringSeconds*state.scoreOps), 10);
				}
			}

			// take the next step
			switch (step) {

				case Energy: {

					// get the next energy conf, if any
					ConfSearch.ScoredConf conf = energyConfs.nextConf();
					if (conf != null) {
						numEnergyConfsEnumerated++;
					}
					if (conf == null || conf.getScore() == Double.POSITIVE_INFINITY) {
						hasEnergyConfs = false;
						keepStepping = false;
						break;
					}

					numConfsEnergied++;

					class EnergyResult {
						ConfSearch.EnergiedConf econf;
						BigDecimal scoreWeight;
						BigDecimal energyWeight;
						Stopwatch stopwatch = new Stopwatch();
					}

					ecalc.tasks.submit(
						() -> {
							// compute one energy and weights (and time it)
							EnergyResult result = new EnergyResult();
							result.stopwatch.start();
							result.econf = ecalc.calcEnergy(conf, confTable);
							result.scoreWeight = bcalc.calc(result.econf.getScore());
							result.energyWeight = bcalc.calc(result.econf.getEnergy());
							result.stopwatch.stop();
							return result;
						},
						(result) -> {
							onEnergy(result.econf, result.scoreWeight, result.energyWeight, result.stopwatch.getTimeS());
						}
					);

					break;
				}

				case Score: {

					// Boolean to get the first score conf and store it
					boolean collectScore = false;
					if (numScoreConfsEnumerated == 0){
						collectScore = true;
					}
					// gather the scores
					List<ConfSearch.ScoredConf> confs = new ArrayList<>();
					for (int i=0; i<numScores; i++) {

						// get the next score conf, if any
						ConfSearch.ScoredConf conf = scoreConfs.nextConf();
						if (conf != null) {
							numScoreConfsEnumerated++;
						}
						if (conf == null || conf.getScore() == Double.POSITIVE_INFINITY) {
							hasScoreConfs = false;
							break;
						}

						confs.add(conf);
					}

					// manually score the first conf to get the first upper bound
					if(collectScore){
						state.firstScoreWeight = bcalc.calc(confs.get(0).getScore());
					}

					class ScoreResult {
						public List<Double> scores = new ArrayList<>();
						List<BigDecimal> scoreWeights = new ArrayList<>();
						Stopwatch stopwatch = new Stopwatch();
					}

					ecalc.tasks.submit(
						() -> {
							// compute the weights (and time it)
							ScoreResult result = new ScoreResult();
							result.stopwatch.start();
							for (ConfSearch.ScoredConf conf : confs) {
								result.scoreWeights.add(bcalc.calc(conf.getScore()));
								result.scores.add(conf.getScore());
							}
							result.stopwatch.stop();
							return result;
						},
						(result) -> {
							onScores(result.scoreWeights, result.stopwatch.getTimeS());
						}
					);

					break;
				}

				case None:
					// out of energy confs and score confs
					// theoretically, this shouldn't happen without hitting our epsilon target, right?
					keepStepping = false;
			}
		}

		// wait for all the scores and energies to come in
		ecalc.tasks.waitForFinish();

		// update the pfunc values from the state
		values.qstar = state.getLowerBound();
		values.qprime = new BigMath(PartitionFunction.decimalPrecision)
			.set(state.getUpperBound())
			.sub(state.getLowerBound())
			.get();

		// we stopped stepping, all the score and energies are accounted for,
		// so update the pfunc status now
		// (allow later status assignments to overwrite previous ones)

		// did we run out of low energies?
		if (!state.hasLowEnergies()) {
			status = Status.OutOfLowEnergies;
		}

		// did we run out of conformations?
		if (!hasEnergyConfs) {
			status = Status.OutOfConformations;
		}

		// did we hit the epsilon target?
		if (state.epsilonReached(targetEpsilon)) {
			status = Status.Estimated;
			System.out.println(String.format("Total Z upper bound reduction through minimizations: %12.6e",state.cumulativeZReduction));
			System.out.println(String.format("Average Z upper bound reduction per minimizations: %12.6e",state.cumulativeZReduction.divide(new BigDecimal(state.numEnergiedConfs),
					new MathContext(BigDecimal.ROUND_HALF_UP))));
		}

		// did we drop below the stability threshold?
		if (!state.isStable(stabilityThreshold)) {
			status = Status.Unstable;
		}
	}

	private void onEnergy(ConfSearch.EnergiedConf econf, BigDecimal scoreWeight, BigDecimal energyWeight, double seconds) {

		synchronized (this) { // don't race the main thread

			// update the state
			state.energyWeightSum = state.energyWeightSum.add(energyWeight);
			state.lowerScoreWeightSum = state.lowerScoreWeightSum.add(scoreWeight);
			state.numEnergiedConfs++;
			state.energyOps = 1.0/seconds;
			if (MathTools.isLessThan(scoreWeight, state.minLowerScoreWeight)) {
				state.minLowerScoreWeight = scoreWeight;
			}

			// set the slope for the energy axis
			double delta = state.calcDelta();
			state.dEnergy = calcSlope(delta, state.prevDelta, state.dScore);
			state.prevDelta = delta;

			state.cumulativeZReduction = state.cumulativeZReduction.add(scoreWeight.subtract(energyWeight));
			int minimizationSize = econf.getAssignments().length;
			if(state.minList.size() < minimizationSize) {
				state.minList.addAll(new ArrayList<Integer>(Collections.nCopies(minimizationSize - state.minList.size(), 0)));
			}
            state.minList.set(minimizationSize-1, state.minList.get(minimizationSize-1)+1);

			// the other direction could be different now, let's be more likely to explore it
			state.dScore *= 2.0;


			// report progress if needed
			if (isReportingProgress) {
				int[] x = econf.getAssignments();
				System.out.println("["+SimpleConfSpace.formatConfRCs(econf)+"] "+String.format("conf:%4d, score:%12.6f, energy:%12.6f, bounds:[%12e,%12e], delta:%.6f, time:%10s, heapMem:%s, extMem:%s",
					state.numEnergiedConfs,
					econf.getScore(), econf.getEnergy(),
					state.getLowerBound().doubleValue(), state.getUpperBound().doubleValue(),
					state.calcDelta(),
					stopwatch.getTime(2),
					JvmMem.getOldPool(),
					ExternalMemory.getUsageReport()
				));
			}

			// update the trace if needed
			if (trace != null) {
				trace.step(state.numScoredConfs, state.numEnergiedConfs, state.calcDelta());
			}
		}

		// report confs if needed
		if (confListener != null) {
			confListener.onConf(econf);
		}
	}

	private void onScores(List<BigDecimal> scoreWeights, double seconds) {

		synchronized (this) { // don't race the main thread
		    // If this is the first score,

			// update the state
			for (BigDecimal weight : scoreWeights) {
				state.upperScoreWeightSum = state.upperScoreWeightSum.add(weight);
				if (MathTools.isLessThan(weight, state.minUpperScoreWeight)) {
					state.minUpperScoreWeight = weight;
				}
			}
			state.numScoredConfs += scoreWeights.size();
			state.scoreOps = scoreWeights.size()/seconds;

			// set the slope for the score axis
			double delta = state.calcDelta();
			state.dScore = calcSlope(delta, state.prevDelta, state.dEnergy);
			state.prevDelta = delta;

			// the other direction could be different now, let's be more likely to explore it
			state.dEnergy *= 2.0;

			// update the trace if needed
			if (trace != null) {
				trace.step(state.numScoredConfs, state.numEnergiedConfs, state.calcDelta());
			}
		}
	}

	private static double calcSlope(double delta, double prevDelta, double otherSlope) {

		double slope = delta - prevDelta;
		// NOTE: this should be 0 or less

		// is the slope totally flat?
		if (slope >= 0.0) {

			// we can't descend on flat slopes, so just pretend it's not flat
			// but we don't want to explore this axis again soon,
			// so set the slope to something a bit flatter than the other axis
			slope = otherSlope/10.0;
		}

		return slope;
	}
	@Override
	public PartitionFunction.Result makeResult() {
	    //Record original bounds
		BigDecimal startLowerBound = BigDecimal.ZERO;
		BigDecimal startUpperBound = state.numConfs.multiply(state.firstScoreWeight);
	    //Record Z reductions
		BigDecimal lowerFullMin = state.getLowerBound(); //Pfunc lower bound improvement from full minimization
		BigDecimal lowerConfUpperBound = BigDecimal.ZERO; //Pfunc lower bound improvement from conf upper bounds, K* has none
		BigDecimal upperFullMin = state.cumulativeZReduction; //Pfunc upper bound improvement from full minimization
		BigDecimal upperPartialMin = BigDecimal.ZERO; //Pfunc upper bound improvement from partial minimization corrections, K* has none

		// first need to calculate upper bound without energied confs
		BigDecimal finalUpperBoundNoEnergies = state.getUpperBoundNoE();
		BigDecimal upperConfLowerBound = startUpperBound.subtract(finalUpperBoundNoEnergies);

		PartitionFunction.Result result = new PartitionFunction.Result(getStatus(), getValues(), getNumConfsEvaluated());
		return result;
	}
}
