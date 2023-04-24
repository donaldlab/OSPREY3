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
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.*;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static edu.duke.cs.osprey.tools.Log.log;


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
public class NewGradientDescentPfunc implements PartitionFunction.WithConfDB, PartitionFunction.WithExternalMemory {

	private static BigMath bigMath() {
		return new BigMath(PartitionFunction.decimalPrecision);
	}

	private static class State {

		BigDecimal numConfs;

		// upper bound (score axis) vars
		long numScoredConfs = 0;
		BigDecimal upperScoreWeightSum = BigDecimal.ZERO;
		BigDecimal minUpperScoreWeight = MathTools.BigPositiveInfinity;

		// lower bound (energy axis) vars
		long numEnergiedConfs = 0;
		BigDecimal lowerScoreWeightSum = BigDecimal.ZERO;
		BigDecimal energyWeightSum = BigDecimal.ZERO;
		BigDecimal minLowerScoreWeight = MathTools.BigPositiveInfinity;
		BigDecimal cumulativeZReduction = BigDecimal.ZERO;
		ArrayList<Integer> minList = new ArrayList<>();
		BigDecimal firstScoreWeight = BigDecimal.ZERO;

		// estimate of inital rates
		// (values here aren't super imporant since they get tuned during execution,
		// but make scores much faster than energies so we don't get a slow start)
		double scoreOps = 100.0;
		double energyOps = 1.0;

		double prevDelta = 1.0;
		double dEnergy = -1.0;
		double dScore = -1.0;

		long lastReportNs = 0;

		State(BigInteger numConfs) {
			this.numConfs = new BigDecimal(numConfs);
		}

		double calcDelta() {
			BigDecimal upperBound = getUpperBound();
			if (MathTools.isZero(upperBound) || MathTools.isInf(upperBound)) {
				return 1.0;
			}
			return bigMath()
					.set(upperBound)
					.sub(getLowerBound())
					.div(upperBound)
					.get()
					.doubleValue();
		}

		public BigDecimal getLowerBound() {
			return energyWeightSum;
		}

		@SuppressWarnings("unused")
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

			return bigMath()

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

			return bigMath()

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
			return String.format("upper: count %d  sum %s  min %s     lower: count %d  score sum %s  energy sum %s",
					numScoredConfs, Log.formatBig(upperScoreWeightSum), Log.formatBig(minUpperScoreWeight),
					numEnergiedConfs, Log.formatBig(lowerScoreWeightSum), Log.formatBig(energyWeightSum)
			);
		}
	}

	private static enum Step {
		None,
		Score,
		Energy
	}


    public final ConfEnergyCalculator ecalc;
	public final BigInteger numConfsBeforePruning;

	private double targetEpsilon = Double.NaN;
	private BigDecimal stabilityThreshold = BigDecimal.ZERO;
	private ConfListener confListener = null;
	private boolean isReportingProgress = false;
	private Stopwatch stopwatch = new Stopwatch().start();
	private ConfSearch scoreConfs;
	private ConfSearch energyConfs;

	private static BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
	private boolean usePreciseBcalc = true;

	private Status status = null;
	private Values values = null;
	private State state = null;

	private boolean hasEnergyConfs = true;
	private boolean hasScoreConfs = true;
	private long numEnergyConfsEnumerated = 0;
	private long numScoreConfsEnumerated = 0;

	private ConfDB confDB = null;
	private ConfDB.Key confDBKey = null;

	private boolean useExternalMemory = false;
	private RCs rcs = null;

	private PfuncSurface surf = null;
	private PfuncSurface.Trace trace = null;

	private final TaskExecutor te;
	private PosInterGen posInterGen;

	public NewGradientDescentPfunc(ConfEnergyCalculator ecalc, ConfSearch upperBoundConfs, ConfSearch lowerBoundConfs, BigInteger numConfsBeforePruning, PosInterGen posInterGen, TaskExecutor te) {
		this.ecalc = ecalc;
		this.scoreConfs = upperBoundConfs;
		this.energyConfs = lowerBoundConfs;
		this.numConfsBeforePruning = numConfsBeforePruning;
		this.posInterGen = posInterGen;
		this.te = te;
	}

	private Integer instanceId = null;

	@Override
	public void setInstanceId(int instanceId) {
		this.instanceId = instanceId;
	}

	public int instanceIdOrThrow() {
		if (instanceId == null) {
			throw new IllegalStateException("no instance ID set, task doesn't know what context to use");
		}
		return instanceId;
	}

	/**
	 * If true, uses BoltzmannCalculator.calcPrecise() instead of BoltzmannCalculator.calc().
	 * The "precise" version of the function has more precision,
	 * and much much better performance on inputs with higher magnitude (eg, -100, -1000).
	 * Generally, this results in significantly faster performance for pfunc calculations.
	 */
	public NewGradientDescentPfunc setPreciseBcalc(boolean val) {
		usePreciseBcalc = val;
		return this;
	}

	ConfDB.ConfTable getConfTable() {
		if (confDB == null) {
			return null;
		}

		return confDB.get(confDBKey);
	}

	BigDecimal bcalc(double val) {
		if (usePreciseBcalc) {
			return bcalc.calcPrecise(val);
		}

		return bcalc.calc(val);
	}

	private void saveToConfDb(ConfSearch.EnergiedConf econf) {
		var confTable = getConfTable();
		if (confTable != null) {
			confTable.setBounds(econf, TimeTools.getTimestampNs());
		}
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

	@Override
	public int getParallelism() {
		return 1;
//		return ecalc.tasks.getParallelism(); // ALWAYS 1
	}

	@Override
	public void setConfDB(ConfDB confDB, ConfDB.Key key) {
		this.confDB = confDB;
		this.confDBKey = key;
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
	public void init(double targetEpsilon) {

		if (targetEpsilon < 0.0) {
			throw new IllegalArgumentException("target epsilon must at least zero");
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

		// split the confs between the upper and lower bounds if needed
		if (energyConfs == null) {
			ConfSearch.Splitter confsSplitter = new ConfSearch.Splitter(scoreConfs, useExternalMemory, rcs);
			scoreConfs = confsSplitter.first;
			energyConfs = confsSplitter.second;
		}
	}

	@Override
	public void setStabilityThreshold(BigDecimal val) {
		this.stabilityThreshold = val;
	}

	static class EnergyResult {
		public ConfSearch.EnergiedConf econf;
		public BigDecimal scoreWeight;
		public BigDecimal energyWeight;
		public Stopwatch stopwatch = new Stopwatch();
	}

	static class ScoreResult {
		public List<Double> scores = new ArrayList<>();
		public List<BigDecimal> scoreWeights = new ArrayList<>();
		public Stopwatch stopwatch = new Stopwatch();
	}

	private EnergyResult computeEnergyResult(ConfSearch.ScoredConf conf) {
		EnergyResult result = new EnergyResult();
		result.stopwatch.start();

		// check if it's already been computed:

		var confSpace = ecalc.confSpace();
		var assignments = conf.getAssignments();

		/*
		var putative = getConfTable().getEnergied(assignments);
		if (putative != null) {
			result.econf = putative;
			result.scoreWeight = bcalc(putative.getScore());
			result.energyWeight = bcalc(putative.getEnergy());

			return result;
		}
		 */

		var inters = posInterGen.all(confSpace, assignments);
		result.econf = ecalc.minimizeEnergy(conf, inters);
		result.scoreWeight = bcalc(result.econf.getScore());
		result.energyWeight = bcalc(result.econf.getEnergy());
		saveToConfDb(result.econf);
		result.stopwatch.stop();

		return result;
	}

	private ScoreResult computeScoreConfs(Iterable<ConfSearch.ScoredConf> confs) {
		// compute the weights (and time it)
		ScoreResult result = new ScoreResult();
		result.stopwatch.start();
		for (ConfSearch.ScoredConf conf : confs) {
			result.scoreWeights.add(bcalc(conf.getScore()));
			result.scores.add(conf.getScore());
		}
		result.stopwatch.stop();
		return result;
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
					te.submit(() -> computeEnergyResult(conf),
							(result) -> onEnergy(result.econf, result.scoreWeight, result.energyWeight, result.stopwatch.getTimeS())
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
							/*
							if (numScoreConfsEnumerated > 10)
							{
								hasScoreConfs = false;
								break;
							}
							 */
						}
						if (conf == null || conf.getScore() == Double.POSITIVE_INFINITY) {
							hasScoreConfs = false;
							break;
						}

						confs.add(conf);
					}

					if (!confs.isEmpty()) {
						te.submit(
							() -> computeScoreConfs(confs),
							(result) -> onScores(result.scoreWeights, result.stopwatch.getTimeS())
						);
					}

					break;
				}

				case None:
					// out of energy confs and score confs
					// theoretically, this shouldn't happen without hitting our epsilon target, right?
					keepStepping = false;
			}
		}

		// wait for all the scores and energies to come in
		te.waitForFinish();

		// update the pfunc values from the state
		values.qstar = state.getLowerBound();
		values.qprime = bigMath()
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
			if (isReportingProgress) {
				log("Total Z upper bound reduction through minimizations: %12.6e", state.cumulativeZReduction);
				log("Average Z upper bound reduction per minimizations: %12.6e", bigMath().set(state.cumulativeZReduction).div(state.numEnergiedConfs).get());
			}
		}

		// did we drop below the stability threshold?
		if (!state.isStable(stabilityThreshold)) {
			status = Status.Unstable;
		}
	}

	private void onEnergy(ConfSearch.EnergiedConf econf, BigDecimal scoreWeight, BigDecimal energyWeight, double seconds) {

		synchronized (this) { // don't race the main thread

			// update the state
			state.energyWeightSum = bigMath()
					.set(state.energyWeightSum)
					.add(energyWeight)
					.get();
			state.lowerScoreWeightSum = bigMath()
					.set(state.lowerScoreWeightSum)
					.add(scoreWeight)
					.get();
			state.numEnergiedConfs++;
			state.energyOps = 1.0/seconds;
			if (MathTools.isLessThan(scoreWeight, state.minLowerScoreWeight)) {
				state.minLowerScoreWeight = scoreWeight;
			}

			// set the slope for the energy axis
			double delta = state.calcDelta();
			state.dEnergy = calcSlope(delta, state.prevDelta, state.dScore);
			state.prevDelta = delta;

			state.cumulativeZReduction = bigMath()
					.set(state.cumulativeZReduction)
					.add(scoreWeight)
					.sub(energyWeight)
					.get();
			int minimizationSize = econf.getAssignments().length;
			if (state.minList.size() < minimizationSize) {
				state.minList.addAll(new ArrayList<>(Collections.nCopies(minimizationSize - state.minList.size(), 0)));
			}
			if (minimizationSize > 0) {
				state.minList.set(minimizationSize-1, state.minList.get(minimizationSize-1)+1);
			}

			// the other direction could be different now, let's be more likely to explore it
			state.dScore *= 2.0;

			// report progress if needed
			if (isReportingProgress) {
				log("[%s] scores:%8d, confs:%4d, score:%12.6f, energy:%12.6f, bounds:[%12f,%12f] (log10p1), delta:%.6f, time:%10s, heapMem:%s",
						SimpleConfSpace.formatConfRCs(econf),
						state.numScoredConfs,
						state.numEnergiedConfs,
						econf.getScore(), econf.getEnergy(),
						MathTools.log10p1(state.getLowerBound()), MathTools.log10p1(state.getUpperBound()),
						state.calcDelta(),
						stopwatch.getTime(2),
						JvmMem.getOldPool()
				);
				state.lastReportNs = System.nanoTime();
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

			// If this is the first score, save it to get the first upper bound
			if (state.numScoredConfs == 0) {
				state.firstScoreWeight = scoreWeights.get(0);
			}

			// update the state
			for (BigDecimal weight : scoreWeights) {
				state.upperScoreWeightSum = bigMath()
						.set(state.upperScoreWeightSum)
						.add(weight)
						.get();
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

			// report progress if needed (but not more than once every second)
			if (isReportingProgress) {
				long nowNs = System.nanoTime();
				if (nowNs - state.lastReportNs > 1_000_000_000L) {
					log("[%s] scores:%8d, confs:%4d, score:%12s, energy:%12s, bounds:[%12f,%12f] (log10p1), delta:%.6f, time:%10s, heapMem:%s, energyOps:%.6f, scoreOps:%.6f",
							String.format("%" + (ecalc.confSpace().numPos()*6 - 1) + "s", ""),
							state.numScoredConfs,
							state.numEnergiedConfs,
							"", "",
							MathTools.log10p1(state.getLowerBound()), MathTools.log10p1(state.getUpperBound()),
							state.calcDelta(),
							stopwatch.getTime(2),
							JvmMem.getOldPool(),
							state.energyOps,
							state.scoreOps
					);
					state.lastReportNs = nowNs;
				}
			}

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
	public Result makeResult() {

		// TODO: all of this seems to have no effect on the result ... should we remove it?
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

		return new Result(getStatus(), getValues(), getNumConfsEvaluated());
	}
}
