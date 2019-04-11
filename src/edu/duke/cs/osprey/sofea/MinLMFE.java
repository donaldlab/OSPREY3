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

package edu.duke.cs.osprey.sofea;


import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.tools.HashCalculator;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;
import edu.duke.cs.osprey.tools.UnpossibleError;
import edu.duke.cs.osprey.tools.resultdoc.Plot;
import edu.duke.cs.osprey.tools.resultdoc.ResultDoc;

import java.io.File;
import java.math.BigInteger;
import java.util.*;
import java.util.function.Function;

import static edu.duke.cs.osprey.tools.Log.log;


/**
 * Finishes SOFEA computation when we've found the best K sequences
 * that minimize an LMFE objective function.
 */
public class MinLMFE implements Sofea.Criterion {

	/**
	 * The multi-state objective function.
	 */
	public final MultiStateConfSpace.LMFE objective;

	/**
	 * The number of best sequences to find.
	 */
	public final int numSequences;

	/**
	 * The best desired precision for a state free energy.
	 */
	public final double minFreeEnergyWidth;

	private final HashSet<StateSeq> finishedSequenced = new HashSet<>();
	private final HashSet<Integer> finishedUnsequenced = new HashSet<>();

	private static class StateSeq {

		final int stateIndex;
		final int[] seq;

		StateSeq(MultiStateConfSpace.State state, Sequence seq) {
			this.stateIndex = state.sequencedIndex;
			this.seq = seq.rtIndices;
		}

		@Override
		public int hashCode() {
			return HashCalculator.combineHashes(
				stateIndex + 1,
				Arrays.hashCode(seq)
			);
		}

		@Override
		public boolean equals(Object obj) {
			return obj instanceof StateSeq && equals((StateSeq)obj);
		}

		public boolean equals(StateSeq other) {
			return this.stateIndex == other.stateIndex
				&& Arrays.equals(this.seq, other.seq);
		}
	}

	public MinLMFE(MultiStateConfSpace.LMFE objective, int numSequences, double minFreeEnergyWidth) {

		// make sure we have enough sequences
		BigInteger maxNumSequences = objective.confSpace.seqSpace.getNumSequences();
		if (BigInteger.valueOf(numSequences).compareTo(maxNumSequences) > 0) {
			throw new IllegalArgumentException(String.format("conf space only has %d sequences, can't find %d",
				maxNumSequences, numSequences
			));
		}

		this.objective = objective;
		this.numSequences = numSequences;
		this.minFreeEnergyWidth = minFreeEnergyWidth;
	}

	public class TopSequences {

		public final List<Sofea.SeqResult> sequences = new ArrayList<>();
		public Sofea.SeqResult nextLowest = null;

		public boolean isTop(DoubleBounds bounds) {
			return nextLowest == null || bounds.upper < nextLowest.lmfeFreeEnergy.upper;
		}

		public void add(Sofea.SeqResult result) {

			assert (isTop(result.lmfeFreeEnergy));

			sequences.add(result);

			// sort the sequences by objective upper bound
			sequences.sort(Comparator.comparing(r -> r.lmfeFreeEnergy.upper));

			// trim down to the best K sequences (allowing ties, but kick out +inf), and update the cutoff
			int numTopK = 0;
			Double lastVal = null;
			for (int i=0; i<sequences.size(); i++) {

				Sofea.SeqResult r = sequences.get(i);
				double val = r.lmfeFreeEnergy.upper;

				// immediately trim +inf out of top K
				if (val == Double.POSITIVE_INFINITY) {
					trimAt(i);
					break;
				}

				if (lastVal == null) {

					// no top K yet results, accept the first result
					numTopK++;
					lastVal = val;

				} else if (val == lastVal) {

					// same as last result in top K, allow ties
					numTopK++;

				} else if (val > lastVal) {

					// worse than last result in top K
					if (numTopK < numSequences) {

						// have room for more in top K, add this result
						numTopK++;
						lastVal = val;

					} else {

						// no more room in top k, trim the rest of the list
						trimAt(i);
						break;
					}
				} else {
					throw new UnpossibleError();
				}
			}
		}

		private void trimAt(int i) {
			nextLowest = sequences.get(i);
			// but don't keep +inf sequences even at the next lowest
			if (nextLowest.lmfeFreeEnergy.upper == Double.POSITIVE_INFINITY) {
				nextLowest = null;
			}
			sequences.subList(i, sequences.size()).clear();
		}
	}

	@Override
	public Filter filterNode(MultiStateConfSpace.State state, int[] conf, BoltzmannCalculator bcalc) {

		// filter out nodes for sequences we've determined sufficiently already!!
		// we're spending too much time minizing confs that don't help us get to the termination criterion
		boolean isFinished;
		if (state.isSequenced) {

			isFinished = finishedSequenced.contains(new StateSeq(
				state,
				state.confSpace.seqSpace.makeSequence(state.confSpace, conf)
			));

		} else {

			isFinished = finishedUnsequenced.contains(state.unsequencedIndex);
		}

		// if the state/sequence is already sufficiently precise, put its nodes back on the queue for later
		if (isFinished) {
			return Filter.Requeue;
		} else {
			return Filter.Process;
		}
	}

	private void updateFinished(MultiStateConfSpace.State state, DoubleBounds freeEnergyBounds) {
		double width = freeEnergyBounds.size();
		if (width <= minFreeEnergyWidth) {
			finishedUnsequenced.add(state.unsequencedIndex);
		}
	}

	private void updateFinished(MultiStateConfSpace.State state, Sequence seq, DoubleBounds freeEnergyBounds) {
		double width = freeEnergyBounds.size();
		if (width <= minFreeEnergyWidth) {
			finishedSequenced.add(new StateSeq(state, seq));
		}
	}

	public TopSequences getTopSequences(SeqDB seqdb, BoltzmannCalculator bcalc) {

		assert (objective.confSpace == seqdb.confSpace);
		MultiStateConfSpace confSpace = seqdb.confSpace;

		TopSequences topSequences = new TopSequences();

		// get the unsequenced z values
		DoubleBounds[] unsequencedFreeEnergy = new DoubleBounds[confSpace.unsequencedStates.size()];
		for (MultiStateConfSpace.State state : confSpace.unsequencedStates) {
			DoubleBounds freeEnergyBounds = bcalc.freeEnergyPrecise(seqdb.getUnsequencedZSumBounds(state));
			updateFinished(state, freeEnergyBounds);
			unsequencedFreeEnergy[state.unsequencedIndex] = freeEnergyBounds;
		}

		// for each sequence and partial sequence encountered so far...
		for (Map.Entry<Sequence,SeqDB.SeqInfo> entry : seqdb.getSequencedZSumBounds()) {

			Sequence seq = entry.getKey();
			SeqDB.SeqInfo seqInfo = entry.getValue();

			DoubleBounds[] stateFreeEnergies = objective.collectFreeEnergies(state -> {
				if (state.isSequenced) {
					return bcalc.freeEnergyPrecise(seqInfo.zSumBounds[state.sequencedIndex]);
				} else {
					return unsequencedFreeEnergy[state.unsequencedIndex];
				}
			});

			// keep track of sequences that are sufficiently precise
			for (MultiStateConfSpace.State state : confSpace.sequencedStates) {
				updateFinished(state, seq, stateFreeEnergies[state.index]);
			}

			// compute bounds on the objective function
			DoubleBounds objectiveBounds = objective.calc().addAll(stateFreeEnergies).bounds;

			// update the top sequences
			if (topSequences.isTop(objectiveBounds)) {

				// yup, get the rest of the free energies if needed
				for (MultiStateConfSpace.State state : confSpace.states) {
					if (stateFreeEnergies[state.index] == null) {
						if (state.isSequenced) {
							stateFreeEnergies[state.index] = bcalc.freeEnergyPrecise(seqInfo.zSumBounds[state.sequencedIndex]);
						} else {
							stateFreeEnergies[state.index] = unsequencedFreeEnergy[state.unsequencedIndex];
						}
					}
				}

				topSequences.add(new Sofea.SeqResult(
					seq,
					objectiveBounds,
					stateFreeEnergies
				));
			}
		}

		return topSequences;
	}

	@Override
	public Satisfied isSatisfied(SeqDB seqdb, FringeDB fringedbLower, FringeDB fringedbUpper, long pass1Step, long pass2Step, BoltzmannCalculator bcalc) {

		assert (objective.confSpace == seqdb.confSpace);
		MultiStateConfSpace confSpace = seqdb.confSpace;

		TopSequences topSequences = getTopSequences(seqdb, bcalc);

		// report progress
		log("lowest %d/%d sequences by LMFE:", topSequences.sequences.size(), numSequences);
		int i = 0;
		Function<Sofea.SeqResult,String> resultToString = result -> {
			StringBuilder buf = new StringBuilder();
			if (result == null) {
				buf.append("(none yet)");
			} else {
				buf.append(String.format("obj=%s w=%9.4f",
					result.lmfeFreeEnergy.toString(4, 9),
					result.lmfeFreeEnergy.size()
				));
				for (MultiStateConfSpace.State state : confSpace.sequencedStates) {
					DoubleBounds stateFreeEnergy = result.stateFreeEnergies[state.index];
					buf.append(String.format("    %s=%s w=%.4f",
						state.name,
						stateFreeEnergy.toString(4, 9),
						stateFreeEnergy.size()
					));
				}
				buf.append(String.format("    [%s]", result.sequence));
			}
			return buf.toString();
		};
		for (Sofea.SeqResult result : topSequences.sequences) {
			log("\tseq %6d   %s", ++i, resultToString.apply(result));
		}
		log("\tnext lowest: %s", resultToString.apply(topSequences.nextLowest));
		for (MultiStateConfSpace.State state : confSpace.unsequencedStates) {
			DoubleBounds g = bcalc.freeEnergyPrecise(seqdb.getUnsequencedZSumBounds(state));
			log("\t%s=%s w=%.4f", state.name, g.toString(4, 9), g.size());
		}

		// all the top K sequences must be ...
		for (Sofea.SeqResult result : topSequences.sequences) {

			// fully assigned
			if (!result.sequence.isFullyAssigned()) {
				return Satisfied.KeepSweeping;
			}

			// have finite bounds
			if (!Double.isFinite(result.lmfeFreeEnergy.lower) || !Double.isFinite(result.lmfeFreeEnergy.upper)) {
				return Satisfied.KeepSweeping;
			}
		}

		// if all possible sequences are already top, we're done
		if (BigInteger.valueOf(topSequences.sequences.size()).compareTo(confSpace.seqSpace.getNumSequences()) >= 0) {
			log("All sequences found, terminating");
			return Satisfied.Terminate;

		}

		// if we don't have enough sequences, keep going
		if (topSequences.nextLowest == null) {
			return Satisfied.KeepSweeping;
		}

		// if all top sequences are sufficiently precise, we're done (because no more refinement is possible)
		Function<Sofea.SeqResult,Boolean> seqFinished = result ->
			confSpace.states.stream().allMatch(state -> {
				if (state.isSequenced) {
					return finishedSequenced.contains(new StateSeq(state, result.sequence));
				} else {
					return finishedUnsequenced.contains(state.unsequencedIndex);
				}
			});
		boolean allFinished = seqFinished.apply(topSequences.nextLowest)
			&& topSequences.sequences.stream().allMatch(result -> seqFinished.apply(result));
		if (allFinished) {
			log("All sequences sufficiently precise, terminating");
			return Satisfied.Terminate;
		}

		// if all top sequences are distinct from the next lowest, we're done
		boolean allDistinct = topSequences.sequences.stream().allMatch(result ->
			result.lmfeFreeEnergy.upper <= topSequences.nextLowest.lmfeFreeEnergy.lower
		);
		if (allDistinct) {
			log("Found top %d sequences, terminating", topSequences.sequences.size());
			return Satisfied.Terminate;
		}

		// nope, keep going
		return Satisfied.KeepSweeping;
	}

	public void makeResultDoc(SeqDB seqdb, File file) {

		BoltzmannCalculator bcalc = new BoltzmannCalculator(seqdb.mathContext);
		TopSequences topSequences = getTopSequences(seqdb, bcalc);

		try (ResultDoc doc = new ResultDoc(file)) {

			doc.h1("SOFEA Results");
			{
				Plot plot = new Plot();
				plot.ylabels = new ArrayList<>();

				Plot.Intervals intervalsMisc = plot.new IntervalsX();
				intervalsMisc.name = "";
				intervalsMisc.data = new ArrayList<>();

				Plot.Intervals[] intervalsSequenced = new Plot.IntervalsX[seqdb.confSpace.sequencedStates.size()];
				for (MultiStateConfSpace.State state : seqdb.confSpace.sequencedStates) {
					intervalsSequenced[state.sequencedIndex] = plot.new IntervalsX();
					intervalsSequenced[state.sequencedIndex].name = state.name;
					intervalsSequenced[state.sequencedIndex].data = new ArrayList<>();
				}

				Plot.Intervals intervalsObjective = plot.new IntervalsX();
				intervalsObjective.name = "objective LMFE";
				intervalsObjective.data = new ArrayList<>();

				// show unsequenced states
				for (MultiStateConfSpace.State state : seqdb.confSpace.unsequencedStates) {

					BigDecimalBounds z = seqdb.getUnsequencedZSumBounds(state);
					DoubleBounds g = bcalc.freeEnergyPrecise(z);

					plot.ylabels.add(state.name);
					intervalsMisc.data.add(new Plot.Interval(g.lower, g.upper));

					// add padding to the other series, so the data and labels align vertically
					intervalsObjective.data.add(null);
					for (Plot.Intervals intervals : intervalsSequenced) {
						intervals.data.add(null);
					}
				}

				// show sequenced states for the top sequences
				for (Sofea.SeqResult result : topSequences.sequences) {

					plot.ylabels.add(result.sequence.toString());

					for (MultiStateConfSpace.State state : seqdb.confSpace.sequencedStates) {
						DoubleBounds g = result.stateFreeEnergies[state.sequencedIndex];
						intervalsSequenced[state.sequencedIndex].data.add(new Plot.Interval(g.lower, g.upper));
					}

					intervalsObjective.data.add(new Plot.Interval(result.lmfeFreeEnergy.lower, result.lmfeFreeEnergy.upper));

					// add padding to the other series, so the data and labels align vertically
					intervalsMisc.data.add(null);
				}

				// show the next lowest
				plot.ylabels.add("next lowest objective LMFE");
				intervalsMisc.data.add(new Plot.Interval(topSequences.nextLowest.lmfeFreeEnergy.lower, topSequences.nextLowest.lmfeFreeEnergy.upper));

				// add padding to the other series, so the data and labels align vertically
				for (Plot.Intervals intervals : intervalsSequenced) {
					intervals.data.add(null);
				}
				intervalsObjective.data.add(null);

				plot.key = "on tmargin horizontal";
				plot.xlabel = "Free Energy (kcal/mol)";
				plot.xlabelrotate = 30.0;
				plot.width = 960;
				plot.height = 70 + plot.ylabels.size()*40;
				doc.plot(plot);
			}

			doc.println();

			// zoom in on just the objective LMFE results
			doc.h2("Objective LMFE results");
			{
				Plot plot = new Plot();
				plot.ylabels = new ArrayList<>();

				Plot.Intervals intervals = plot.new IntervalsX();
				intervals.name = "objective LMFE";
				intervals.data = new ArrayList<>();

				for (Sofea.SeqResult result : topSequences.sequences) {
					plot.ylabels.add(result.sequence.toString());
					intervals.data.add(new Plot.Interval(result.lmfeFreeEnergy.lower, result.lmfeFreeEnergy.upper));
				}

				// show the next lowest
				plot.ylabels.add("next lowest sequence");
				intervals.data.add(new Plot.Interval(topSequences.nextLowest.lmfeFreeEnergy.lower, topSequences.nextLowest.lmfeFreeEnergy.upper));

				plot.key = "on tmargin horizontal";
				plot.xlabel = "Free Energy (kcal/mol)";
				plot.xlabelrotate = 30.0;
				plot.width = 960;
				plot.height = 70 + plot.ylabels.size()*40;
				doc.plot(plot);
			}
		}
	}
}
