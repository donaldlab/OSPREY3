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

package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.AStarProgress;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.UpperLowerAStarOrder;
import edu.duke.cs.osprey.astar.conf.pruning.AStarPruner;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.PfuncSurface;
import edu.duke.cs.osprey.markstar.MARKStar;
import edu.duke.cs.osprey.markstar.MARKStarProgress;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.JvmMem;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.*;
import java.util.concurrent.PriorityBlockingQueue;


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
public class GradientDescentMARKStarPfunc implements PartitionFunction.WithConfTable, PartitionFunction.WithExternalMemory {

	private final AStarScorer gScorer;
	private final AStarScorer rigidgScorer;
	private final AStarScorer hScorer;
	private final AStarScorer negatedHScorer;
	private final EnergyMatrix minimizingEnergyMatrix;
	private final EnergyMatrix rigidEnergyMatrix;
	private AStarOrder order;
	private AStarPruner pruner;
	private MARKStarProgress progress;

	public void setCorrections(UpdatingEnergyMatrix corrections) {
		this.correctionMatrix= corrections;
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
		ArrayList<Integer> minList = new ArrayList<Integer>();

		// estimate of inital rates
		// (values here aren't super imporant since they get tuned during execution,
		// but make scores much faster than energies so we don't get a slow start)
		double scoreOps = 100.0;
		double energyOps = 1.0;

		double prevDelta = 1.0;
		double dEnergy = -1.0;
		double dScore = -1.0;
		private MARKStarNode root;

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
			return root.getUpperBound().subtract(lowerScoreWeightSum).add(energyWeightSum);
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
	private MARKStarNode root;
	private Queue<MARKStarNode> queue = new PriorityBlockingQueue<>();
	private UpdatingEnergyMatrix correctionMatrix;

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

	public GradientDescentMARKStarPfunc(SimpleConfSpace confSpace, EnergyMatrix rigidEnergyMatrix,
										EnergyMatrix minimizingEnergyMatrix, RCs rcs,
										ConfEnergyCalculator ecalc) {
		this.ecalc = ecalc;
		this.minimizingEnergyMatrix = minimizingEnergyMatrix;
		this.rigidEnergyMatrix = rigidEnergyMatrix;
		MARKStarNode.ScorerFactory gScorerFactory = (emats) -> new PairwiseGScorer(emats);
		MARKStarNode.ScorerFactory hScorerFactory = (emats) -> new TraditionalPairwiseHScorer(minimizingEnergyMatrix, rcs);
		gScorer = new PairwiseGScorer(minimizingEnergyMatrix);
		rigidgScorer = new PairwiseGScorer(rigidEnergyMatrix);
		negatedHScorer = hScorerFactory.make(new NegatedEnergyMatrix(confSpace, rigidEnergyMatrix));
		hScorer = hScorerFactory.make(minimizingEnergyMatrix);
		EnergyMatrix minimizingEmat = minimizingEnergyMatrix;
		EnergyMatrix rigidEmat = rigidEnergyMatrix;
		root = MARKStarNode.makeRoot(confSpace, rigidEmat, minimizingEmat, rcs,
				gScorerFactory.make(minimizingEmat), hScorerFactory.make(minimizingEmat),
				gScorerFactory.make(rigidEmat),
				new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs), true);
        order = new UpperLowerAStarOrder();
        this.rcs = rcs;
        this.correctionMatrix = new UpdatingEnergyMatrix(confSpace, minimizingEnergyMatrix);
		this.progress = new MARKStarProgress(rcs.getNumPos());
	}

	private void boundLowestBoundConfUnderNode(MARKStarNode startNode, List<MARKStarNode> generatedNodes) {
		Comparator<MARKStarNode> confBoundComparator = Comparator.comparingDouble(o -> o.getConfSearchNode().getConfLowerBound());
		PriorityQueue<MARKStarNode> drillQueue = new PriorityQueue<>(confBoundComparator);
		drillQueue.add(startNode);

		List<MARKStarNode> newNodes = new ArrayList<>();
		int numNodes = 0;
		Stopwatch leafLoop = new Stopwatch().start();
		Stopwatch overallLoop = new Stopwatch().start();
		while(!drillQueue.isEmpty()) {
			numNodes++;
			MARKStarNode curNode = drillQueue.poll();
			MARKStarNode.Node node = curNode.getConfSearchNode();
			ConfIndex index = new ConfIndex(rcs.getNumPos());
			node.index(index);

			if (node.getLevel() < rcs.getNumPos()) {
				MARKStarNode nextNode = drillDown(newNodes, curNode, node);
				newNodes.remove(nextNode);
				drillQueue.add(nextNode);
			}
			else {
				newNodes.add(curNode);
				break;
			}

			//debugHeap(drillQueue, true);
			if(leafLoop.getTimeS() > 10) {
				leafLoop.stop();
				leafLoop.reset();
				leafLoop.start();
				System.out.println(String.format("Processed %d, %s so far. Bounds are now [%12.6e,%12.6e]",numNodes, overallLoop.getTime(2),root.getLowerBound(),root.getUpperBound()));
			}
		}
		generatedNodes.addAll(newNodes);

	}
    private MARKStarNode drillDown(List<MARKStarNode> newNodes, MARKStarNode curNode, MARKStarNode.Node node) {
            ConfIndex confIndex = new ConfIndex(rcs.getNumPos());
            node.index(confIndex);
            // which pos to expand next?
            int nextPos = order.getNextPos(confIndex, rcs);
            assert (!confIndex.isDefined(nextPos));
            assert (confIndex.isUndefined(nextPos));

            // score child nodes with tasks (possibly in parallel)
            List<MARKStarNode> children = new ArrayList<>();
            double bestChildLower = Double.POSITIVE_INFINITY;
            MARKStarNode bestChild = null;
            for (int nextRc : rcs.get(nextPos)) {

                if (hasPrunedPair(confIndex, nextPos, nextRc)) {
                    continue;
                }

                // if this child was pruned dynamically, then don't score it
                if (pruner != null && pruner.isPruned(node, nextPos, nextRc)) {
                    continue;
                }
                Stopwatch partialTime = new Stopwatch().start();
                MARKStarNode.Node child = node.assign(nextPos, nextRc);
                double confLowerBound = Double.POSITIVE_INFINITY;

                // score the child node differentially against the parent node
                if (child.getLevel() < rcs.getNumPos()) {
                    //double confCorrection = correctionMatrix.confE(child.assignments);
                    double diff = gScorer.calcDifferential(confIndex, rcs, nextPos, nextRc);//confCorrection;
                    double rigiddiff = rigidgScorer.calcDifferential(confIndex, rcs, nextPos, nextRc);
                    double hdiff = hScorer.calcDifferential(confIndex, rcs, nextPos, nextRc);
                    double maxhdiff = -negatedHScorer.calcDifferential(confIndex, rcs, nextPos, nextRc);
                    child.gscore = diff;
                    //Correct for incorrect gscore.
                    rigiddiff = rigiddiff - node.gscore + node.rigidScore;
                    child.rigidScore = rigiddiff;

                    confLowerBound = child.gscore + hdiff;
                    double confUpperbound = rigiddiff + maxhdiff;
                    child.computeNumConformations(rcs);
                    /*
                    double lowerbound = minimizingEnergyMatrix.confE(child.assignments);
                    if (diff < confCorrection) {
                        recordCorrection(confLowerBound, confCorrection - diff);
                        confLowerBound = confCorrection + hdiff;
                    }
                    */
                    child.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperbound);
                    progress.reportInternalNode(child.level, child.gscore, child.getHScore(), queue.size(), children.size(), state.prevDelta);
                }
                if (child.getLevel() == rcs.getNumPos()) {
                    double confRigid = rigidgScorer.calcDifferential(confIndex, rcs, nextPos, nextRc);
                    confRigid = confRigid - node.gscore + node.rigidScore;

                    child.computeNumConformations(rcs); // Shouldn't this always eval to 1, given that we are looking at leaf nodes?
                    double confCorrection = gScorer.calcDifferential(confIndex, rcs, nextPos, nextRc);// correctionMatrix.confE(child.assignments);
					/*
                    double lowerbound = minimizingEnergyMatrix.confE(child.assignments);
                    if (lowerbound < confCorrection) {
                        recordCorrection(lowerbound, confCorrection - lowerbound);
                    }
                    */
                    child.setBoundsFromConfLowerAndUpper(confCorrection, confRigid);
                    child.gscore = child.getConfLowerBound();
                    //confLowerBound = lowerbound;
                    child.rigidScore = confRigid;
                    state.numScoredConfs++;
                    progress.reportLeafNode(child.gscore, queue.size(), state.prevDelta);
                }
                partialTime.stop();


                if (Double.isNaN(child.rigidScore))
                    System.out.println("Huh!?");
                MARKStarNode MARKStarNodeChild = curNode.makeChild(child);
                MARKStarNodeChild.markUpdated();
                /*
                if (confLowerBound < bestChildLower) {
                    bestChild = MARKStarNodeChild;
                }
                */
                // collect the possible children
                if (MARKStarNodeChild.getConfSearchNode().getConfLowerBound() < 0) {
                    children.add(MARKStarNodeChild);
                }
                newNodes.add(MARKStarNodeChild);

            }
            return bestChild;
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
		energyConfs = confSearch;
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
		state.root = root;
		values = Values.makeFullRange();
		// don't explicitly check the pruned confs, just lump them together with the un-enumerated confs
		values.pstar = BigDecimal.ZERO;

		hasEnergyConfs = true;
		hasScoreConfs = true;
		numEnergyConfsEnumerated = 0;
		numScoreConfsEnumerated = 0;

		List<MARKStarNode> newNodes = new ArrayList<>();
		queue.add(root);
		//boundLowestBoundConfUnderNode(root,newNodes);
		//queue.addAll(newNodes);
	}

	@Override
	public void setStabilityThreshold(BigDecimal val) {
		this.stabilityThreshold = val;
	}

	public void setRCs(RCs rcs){
		this.rcs = rcs;
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
					numScores = Math.max((int)(scoringSeconds*state.scoreOps), 1);
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

					// gather the scores
					List<MARKStarNode> nodes = new ArrayList<>();
					for(int i = 0; i < numScores; i++) {
					    MARKStarNode conf = queue.poll();//getLeaf();
					    if(conf == null)
					    	continue;
					    if(queue.isEmpty())
					    	processNode(conf);
					    else
							nodes.add(conf);

						// get the next score conf, if any
						if (conf != null) {
							numScoreConfsEnumerated++;
						}
						if (conf == null || conf.getConfSearchNode().getConfLowerBound() == Double.POSITIVE_INFINITY) {
							hasScoreConfs = false;
							break;
						}
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
							for (MARKStarNode node: nodes) {
							    processNode(node);
							    node.markUpdated();
							}
							result.stopwatch.stop();
							return result;
						},
						(result) -> {
							onScores(nodes.size(), result.stopwatch.getTimeS());
						}
					);
					//Debug lines. If you pulled from the repo and see this you can delete it.
					if(false) {
						String bounds = state.getLowerBound()+","+state.getUpperBound();
						if(!MathTools.isInf(state.getUpperBound()))
							bounds = String.format("%12e+%12e", state.getLowerBound(), state.getUpperBound().subtract(state.getLowerBound()));
						System.out.println("Score weights are now: [" + bounds + "]");
						state.printBoundStats();
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
			//Debug printline. Delete if you see it.
			if(false) {
				System.out.println(String.format("Partition function approximation complete: [%12.6e,%12.6e]", state.getLowerBound(), state.getUpperBound()));
				System.out.println(String.format("Score breakdown: q* = %12.6e, q' = %12.6e, p' = %12.6e", values.qstar, values.qprime, values.pstar));
			}

		}

		// did we drop below the stability threshold?
		if (!state.isStable(stabilityThreshold)) {
			status = Status.Unstable;
		}
	}

	private void processNode(MARKStarNode node) {
	    if(node.getConfSearchNode().getLevel() == rcs.getNumPos())
	        processFullConfNode(node);
	    else {
	    	List<MARKStarNode> nodes = new ArrayList<>();
	    	processPartialConfNode(node);
	        //boundLowestBoundConfUnderNode(node, nodes);
	        //queue.addAll(nodes);
		}
	}

	private void processFullConfNode(MARKStarNode node) {
	}

	private void processPartialConfNode(MARKStarNode curNode) {
        // which pos to expand next?
            ConfIndex confIndex = new ConfIndex(rcs.getNumPos());
            MARKStarNode.Node node = curNode.getConfSearchNode();
            node.index(confIndex);
            int nextPos = order.getNextPos(confIndex, rcs);
            assert (!confIndex.isDefined(nextPos));
            assert (confIndex.isUndefined(nextPos));

            //preminimizePartialAsync(curNode);

            // score child nodes with tasks (possibly in parallel)
            List<MARKStarNode> children = new ArrayList<>();
            for (int nextRc : rcs.get(nextPos)) {

                if (hasPrunedPair(confIndex, nextPos, nextRc)) {
                    continue;
                }

                // if this child was pruned dynamically, then don't score it
                if (pruner != null && pruner.isPruned(node, nextPos, nextRc)) {
                    continue;
                }


                Stopwatch partialTime = new Stopwatch().start();
                node.index(confIndex);
                MARKStarNode.Node child = node.assign(nextPos, nextRc);

                // score the child node differentially against the parent node
                if (child.getLevel() < rcs.getNumPos()) {
                    double confCorrection = correctionMatrix.confE(child.assignments);
                    double diff = confCorrection;
                    double rigiddiff = rigidgScorer.calcDifferential(confIndex, rcs, nextPos, nextRc);
                    double hdiff = hScorer.calcDifferential(confIndex, rcs, nextPos, nextRc);
                    double maxhdiff = -negatedHScorer.calcDifferential(confIndex, rcs, nextPos, nextRc);
                    child.gscore = diff;
                    //Correct for incorrect gscore.
                    rigiddiff = rigiddiff - node.gscore + node.rigidScore;
                    child.rigidScore = rigiddiff;

                    double confLowerBound = child.gscore + hdiff;
                    double confUpperbound = rigiddiff + maxhdiff;
                    child.computeNumConformations(rcs);
                    double lowerbound = minimizingEnergyMatrix.confE(child.assignments);
                    if (diff < confCorrection) {
                        recordCorrection(confLowerBound, confCorrection - diff);
                        confLowerBound = confCorrection + hdiff;
                    }
                    child.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperbound);
                    progress.reportInternalNode(child.level, child.gscore, child.getHScore(), queue.size(), children.size(), state.prevDelta);
                }
                if (child.getLevel() == rcs.getNumPos()) {
                    double confRigid = rigidgScorer.calcDifferential(confIndex, rcs, nextPos, nextRc);
                    confRigid = confRigid - node.gscore + node.rigidScore;

                    child.computeNumConformations(rcs); // Shouldn't this always eval to 1, given that we are looking at leaf nodes?
                    double confCorrection = correctionMatrix.confE(child.assignments);
                    double lowerbound = minimizingEnergyMatrix.confE(child.assignments);
                    if (lowerbound < confCorrection) {
                        recordCorrection(lowerbound, confCorrection - lowerbound);
                    }
                    child.setBoundsFromConfLowerAndUpper(confCorrection, confRigid);
                    child.gscore = child.getConfLowerBound();
                    child.rigidScore = confRigid;
                    state.numScoredConfs++;
                    progress.reportLeafNode(child.gscore, queue.size(), state.prevDelta);
                }
                partialTime.stop();


                if (Double.isNaN(child.rigidScore))
                    System.out.println("Huh!?");
                MARKStarNode MARKStarNodeChild = curNode.makeChild(child);
                // collect the possible children
                if (MARKStarNodeChild.getConfSearchNode().getConfLowerBound() < 0) {
                    children.add(MARKStarNodeChild);
                }
                if (!child.isMinimized()) {
                    queue.add(MARKStarNodeChild);
                } else
                    MARKStarNodeChild.computeEpsilonErrorBounds();

            }
            curNode.markUpdated();
        }

	private void recordCorrection(double lowerbound, double v) {
	}

	private boolean hasPrunedPair(ConfIndex confIndex, int nextPos, int nextRc) {
	    return false;
	}


	private void onEnergy(ConfSearch.EnergiedConf econf, BigDecimal scoreWeight, BigDecimal energyWeight, double seconds) {

		synchronized (this) { // don't race the main thread


			if(energyWeight.compareTo(scoreWeight) == 1) {
				System.err.println("Bounds are incorrect:" + energyWeight + " > "
						+ scoreWeight);
			}
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

	private void onScores(int numScored, double seconds) {

		synchronized (this) { // don't race the main thread
		    BigDecimal oldUpper = state.getUpperBound();
			root.computeEpsilonErrorBounds();
			BigDecimal newUpper = state.getUpperBound();
			System.out.println(String.format("Upper bound was %12.6e, now is %12.6e",oldUpper, newUpper));
		    double newEpsilon = state.calcDelta();

			// update the state
			state.numScoredConfs += numScored;
			state.scoreOps = numScored/seconds;

			// set the slope for the score axis
			double delta = newEpsilon;
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
	public Result makeResult() {
		return new Result(getStatus(), getValues(), getNumConfsEvaluated());
	}
}
