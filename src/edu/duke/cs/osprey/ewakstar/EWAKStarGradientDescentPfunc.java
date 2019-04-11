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

package edu.duke.cs.osprey.ewakstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PfuncSurface;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.JvmMem;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;


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
public class EWAKStarGradientDescentPfunc implements EWAKStarPartitionFunction.WithConfTable, EWAKStarPartitionFunction.WithExternalMemory {

	private static class State {

		ArrayList<EnergyCalculator.EnergiedParametricMolecule> epMols = new ArrayList<>();
		HashMap<Double, ConfSearch.ScoredConf> sConfs = new HashMap<>();
		BigDecimal numConfs;

		double curGMEC = Double.POSITIVE_INFINITY;

		double curScore = Double.POSITIVE_INFINITY;
		// upper bound (score axis) vars
		long numScoredConfs = 0;
		BigDecimal upperScoreWeightSum = BigDecimal.ZERO;
		BigDecimal minUpperScoreWeight = MathTools.BigPositiveInfinity;

		// lower bound (energy axis) vars
		long numEnergiedConfs = 0;
		BigDecimal lowerScoreWeightSum = BigDecimal.ZERO;
		BigDecimal energyWeightSum = BigDecimal.ZERO;
		BigDecimal minLowerScoreWeight = MathTools.BigPositiveInfinity;

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

		public void setGMECEnergy(double newEnergy){
			this.curGMEC = newEnergy;
		}

		double calcDiff(){
			return (curScore - curGMEC);
		}

		double calcDelta() {
			BigDecimal upperBound = getUpperBound();
			if (MathTools.isZero(upperBound) || MathTools.isInf(upperBound)) {
				return 1.0;
			}
			return new BigMath(EWAKStarPartitionFunction.decimalPrecision)
				.set(upperBound)
				.sub(getLowerBound())
				.div(upperBound)
				.get()
				.doubleValue();
		}

		public BigDecimal getLowerBound() { return energyWeightSum; }

		public BigDecimal getUpperBound() {
			return new BigMath(EWAKStarPartitionFunction.decimalPrecision)

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

		boolean epsilonReached(double targetEpsilon) {return calcDelta() <= targetEpsilon;}

		boolean energyReached(double targetEnergy) {return calcDiff() >= targetEnergy;}

		boolean confMaxReached(int highestNumConfs) {return numEnergiedConfs >= highestNumConfs;}

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
	private double targetEnergy = 0.0;
	private BigDecimal stabilityThreshold = BigDecimal.ZERO;
	private ConfListener confListener = null;
	private boolean isReportingProgress = false;
	private Stopwatch stopwatch = new Stopwatch().start();
	private ConfSearch scoreConfs = null;
	private ConfSearch energyConfs = null;
	private BoltzmannCalculator bcalc = new BoltzmannCalculator(EWAKStarPartitionFunction.decimalPrecision);

	private Status status = null;
	private Values values = null;
	private State state = null;

	private boolean hasEnergyConfs = true;
	private boolean hasScoreConfs = true;
	private long numEnergyConfsEnumerated = 0;
	private long numScoreConfsEnumerated = 0;
	private int highestNumConfs = 5000;

	private boolean printPDBs = false;

	private ConfDB.ConfTable confTable = null;

	private boolean useExternalMemory = false;
	private RCs rcs = null;

	private PfuncSurface surf = null;
	private PfuncSurface.Trace trace = null;

	public EWAKStarGradientDescentPfunc(ConfEnergyCalculator ecalc) {this.ecalc = ecalc;}

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
	public ArrayList<EnergyCalculator.EnergiedParametricMolecule> getEpMols(){
		return state.epMols;
	}

	@Override
	public HashMap<Double, ConfSearch.ScoredConf> getSConfs(){
		return state.sConfs;
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
	public void init(ConfSearch scoreConfs, ConfSearch energyConfs, BigInteger numConfsBeforePruning, double targetEpsilon, double targetEnergy, int highestNumConfs, boolean printPDBs) {

		if (targetEpsilon <= 0.0 || targetEnergy < 0) {
			throw new IllegalArgumentException("target epsilon and target energy must be greater than zero");
		}

		this.printPDBs = printPDBs;
		this.targetEpsilon = targetEpsilon;
		this.targetEnergy = targetEnergy;

		this.highestNumConfs = highestNumConfs;
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

		this.scoreConfs = scoreConfs;
		this.energyConfs = energyConfs;
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
					&& state.hasLowEnergies()
					&& !state.energyReached(targetEnergy)
					&& !state.confMaxReached(highestNumConfs);
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
						EnergyCalculator.EnergiedParametricMolecule epmol;
						BigDecimal scoreWeight;
						BigDecimal energyWeight;
						Stopwatch stopwatch = new Stopwatch();
					}

					ecalc.tasks.submit(
						() -> {
							// compute one energy and weights (and time it)
							EnergyResult result = new EnergyResult();
							result.stopwatch.start();
							if(printPDBs){result.epmol = ecalc.calcEnergy(new RCTuple(conf.getAssignments()));}
							result.econf = ecalc.calcEnergy(conf, confTable);
							if (printPDBs) { state.sConfs.put(result.econf.getEnergy(), conf);}
							result.scoreWeight = bcalc.calc(result.econf.getScore());
							result.energyWeight = bcalc.calc(result.econf.getEnergy());
							result.stopwatch.stop();
							return result;
						},
						(result) -> {
							onEnergy(result.epmol, result.econf, result.scoreWeight, result.energyWeight, result.stopwatch.getTimeS());
						}
					);

					break;
				}

				case Score: {

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

					class ScoreResult {
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
		values.qprime = new BigMath(EWAKStarPartitionFunction.decimalPrecision)
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
			status = Status.EpsilonReached;
		}

		//did we hit the energy window?
		if(state.energyReached(targetEnergy)){
			status = Status.EnergyReached;
		}

		//did we hit the conformation limit?
		if(state.numEnergiedConfs >= highestNumConfs){
			status = Status.ConfLimitReached;
		}

		// did we drop below the stability threshold?
		if (!state.isStable(stabilityThreshold)) {
			status = Status.Unstable;
		}
	}

	private void onEnergy(EnergyCalculator.EnergiedParametricMolecule epmol, ConfSearch.EnergiedConf econf, BigDecimal scoreWeight, BigDecimal energyWeight, double seconds) {

		synchronized (this) { // don't race the main thread

			// update the state
			state.energyWeightSum = state.energyWeightSum.add(energyWeight);
			state.epMols.add(epmol);

			if (state.curGMEC >= econf.getEnergy())
				state.setGMECEnergy(econf.getEnergy());

			state.curScore = econf.getScore();

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

			// the other direction could be different now, let's be more likely to explore it
			state.dScore *= 2.0;

			// report progress if needed
			if (isReportingProgress) {
				System.out.println(String.format("conf:%4d, score:%12.6f, energy:%12.6f, bounds:[%12e,%12e], delta:%.6f, time:%10s, heapMem:%s, extMem:%s",
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
}
