package edu.duke.cs.osprey.sofea;


import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;
import edu.duke.cs.osprey.tools.UnpossibleError;
import edu.duke.cs.osprey.tools.resultdoc.Plot;
import edu.duke.cs.osprey.tools.resultdoc.ResultDoc;

import java.io.File;
import java.math.BigInteger;
import java.util.*;

import static edu.duke.cs.osprey.tools.Log.log;
import static edu.duke.cs.osprey.tools.Log.logf;


/**
 * Finishes SOFEA computation when we've found the lowest K sequences by LMFE
 */
public class MinLMFE implements Sofea.Criterion {

	public final MultiStateConfSpace.LMFE objective;
	public final int numSequences;
	public final BoltzmannCalculator bcalc;

	public MinLMFE(MultiStateConfSpace.LMFE objective, int numSequences, BoltzmannCalculator bcalc) {

		// make sure we have enough sequences
		BigInteger maxNumSequences = objective.confSpace.seqSpace.getNumSequences();
		if (BigInteger.valueOf(numSequences).compareTo(maxNumSequences) > 0) {
			throw new IllegalArgumentException(String.format("conf space only has %d sequences, can't find %d",
				maxNumSequences, numSequences
			));
		}

		this.objective = objective;
		this.numSequences = numSequences;
		this.bcalc = bcalc;
	}

	public class TopSequences {

		public final List<Sofea.SeqResult> sequences = new ArrayList<>();
		public Double nextLowest = null;

		public boolean isTop(DoubleBounds bounds) {
			return nextLowest == null || bounds.lower < nextLowest;
		}

		public void add(Sofea.SeqResult result) {

			assert (isTop(result.lmfeFreeEnergy));

			sequences.add(result);

			// sort the sequences by objective lower bound
			sequences.sort(Comparator.comparing(r -> r.lmfeFreeEnergy.lower));

			// trim down to the best K sequences (allowing ties, but kick out +inf), and update the cutoff
			int numTopK = 0;
			Double lastVal = null;
			for (int i=0; i<sequences.size(); i++) {

				Sofea.SeqResult r = sequences.get(i);
				double val = r.lmfeFreeEnergy.lower;

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
			nextLowest = sequences.get(i).lmfeFreeEnergy.lower;
			sequences.subList(i, sequences.size()).clear();
		}
	}

	@Override
	public Filter filterNode(MultiStateConfSpace.State state, int[] conf, BoltzmannCalculator bcalc) {
		// TODO: filter out nodes for sequences we've determined sufficiently already?
		// accept all by default
		return Filter.Process;
	}

	public TopSequences getTopSequences(SeqDB seqdb) {

		assert (objective.confSpace == seqdb.confSpace);
		MultiStateConfSpace confSpace = seqdb.confSpace;

		TopSequences topSequences = new TopSequences();

		// get the unsequenced z values
		DoubleBounds[] unsequencedFreeEnergy = new DoubleBounds[confSpace.unsequencedStates.size()];
		for (MultiStateConfSpace.State state : confSpace.unsequencedStates) {
			unsequencedFreeEnergy[state.unsequencedIndex] = bcalc.freeEnergyPrecise(seqdb.getUnsequencedZSumBounds(state));
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
	public Satisfied isSatisfied(SeqDB seqdb, FringeDB fringedb, long sweepCount, BoltzmannCalculator bcalc) {

		assert (objective.confSpace == seqdb.confSpace);
		MultiStateConfSpace confSpace = seqdb.confSpace;

		TopSequences topSequences = getTopSequences(seqdb);

		// report progress
		log("lowest %d/%d sequences by LMFE:", topSequences.sequences.size(), numSequences);
		int i = 0;
		for (Sofea.SeqResult result : topSequences.sequences) {
			logf("\tseq %2d  obj=%s w=%9.4f",
				++i,
				result.lmfeFreeEnergy.toString(4, 9),
				result.lmfeFreeEnergy.size()
			);
			for (MultiStateConfSpace.State state : confSpace.sequencedStates) {
				DoubleBounds stateFreeEnergy = result.stateFreeEnergies[state.index];
				logf("    %s=%s w=%.4f",
					state.name,
					stateFreeEnergy.toString(4, 9),
					stateFreeEnergy.size()
				);
			}
			log("    [%s]", result.sequence);
		}
		log("\tnext lowest: %9.4f", topSequences.nextLowest);
		for (MultiStateConfSpace.State state : confSpace.unsequencedStates) {
			DoubleBounds g = bcalc.freeEnergyPrecise(seqdb.getUnsequencedZSumBounds(state));
			log("\t%s=%s w=%.4f", state.name, g.toString(4, 9), g.size());
		}

		// all the top K sequences must be ...
		for (Sofea.SeqResult result : topSequences.sequences) {

			// fully assigned
			if (!result.sequence.isFullyAssigned()) {
				return Satisfied.KeepIterating;
			}

			// have finite bounds
			if (!Double.isFinite(result.lmfeFreeEnergy.lower) || !Double.isFinite(result.lmfeFreeEnergy.upper)) {
				return Satisfied.KeepIterating;
			}
		}

		// if all possible sequences are already top, we're done
		if (BigInteger.valueOf(topSequences.sequences.size()).compareTo(confSpace.seqSpace.getNumSequences()) >= 0) {
			return Satisfied.Terminate;

		// if we don't have enough sequences, keep going
		} else if (topSequences.nextLowest == null) {
			return Satisfied.KeepIterating;
		}

		// if any top sequence overlaps the next lowest, keep going
		for (Sofea.SeqResult result : topSequences.sequences) {
			if (result.lmfeFreeEnergy.upper >= topSequences.nextLowest) {
				return Satisfied.KeepIterating;
			}
		}

		// we're all done here!
		return Satisfied.Terminate;
	}

	public void makeResultDoc(SeqDB seqdb, File file) {

		TopSequences topSequences = getTopSequences(seqdb);

		BoltzmannCalculator bcalc = new BoltzmannCalculator(seqdb.mathContext);
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
				intervalsMisc.data.add(new Plot.Interval(topSequences.nextLowest, topSequences.nextLowest));

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
				intervals.data.add(new Plot.Interval(topSequences.nextLowest, topSequences.nextLowest));

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
