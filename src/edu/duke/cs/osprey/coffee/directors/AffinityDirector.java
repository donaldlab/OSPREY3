package edu.duke.cs.osprey.coffee.directors;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.Coffee;
import edu.duke.cs.osprey.coffee.FreeEnergyCalculator;
import edu.duke.cs.osprey.coffee.NodeProcessor;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.coffee.seqdb.SeqFreeEnergies;
import edu.duke.cs.osprey.coffee.seqdb.SeqInfo;
import edu.duke.cs.osprey.coffee.seqdb.StateZ;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.parallelism.ThreadTools;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;


/**
 * Directs COFFEE to find the K best sequences by binding affinity
 *
 * The complex and design should be mutable, the target should not.
 */
public class AffinityDirector implements Coffee.Director {

	public static class Builder {

		private final MultiStateConfSpace confSpace;
		private final MultiStateConfSpace.State complex;
		private final MultiStateConfSpace.State design;
		private final MultiStateConfSpace.State target;

		private int K = 5;
		private double gWidthMax = 1.0;
		private Timing timing = Timing.Efficient;

		public Builder(MultiStateConfSpace confSpace, MultiStateConfSpace.State complex, MultiStateConfSpace.State design, MultiStateConfSpace.State target) {

			// complex and design states should be sequenced, target state should not
			if (!complex.isSequenced) {
				throw new IllegalArgumentException("The complex state should be sequenced");
			}
			if (!design.isSequenced) {
				throw new IllegalArgumentException("The design state should be sequenced");
			}
			if (target.isSequenced) {
				throw new IllegalArgumentException("The target state should not be sequenced");
			}

			this.confSpace = confSpace;
			this.complex = complex;
			this.design = design;
			this.target = target;
		}

		public Builder(MultiStateConfSpace confSpace, String complexName, String designName, String targetName) {
			this(confSpace, confSpace.getState(complexName), confSpace.getState(designName), confSpace.getState(targetName));
		}

		public Builder setK(int val) {
			K = val;
			return this;
		}

		/**
		 * Sets the largest precision desired for free energy calculations.
		 * Resulting free energy values may be more precise than this maximum value.
		 * If the computation is stopped early, the free energy values may be less precise than this value.
		 */
		public Builder setGWidthMax(double val) {
			gWidthMax = val;
			return this;
		}

		public Builder setTiming(Timing val) {
			timing = val;
			return this;
		}

		public AffinityDirector build() {
			return new AffinityDirector(
				confSpace, complex, design, target,
				K, gWidthMax, timing
			);
		}
	}

	public final MultiStateConfSpace confSpace;
	public final MultiStateConfSpace.State complex;
	public final MultiStateConfSpace.State design;
	public final MultiStateConfSpace.State target;
	public final int K;
	public final double gWidthMax;
	public final Timing timing;

	public final List<SeqFreeEnergies> bestSeqs;

	private final BestSet<SeqFreeEnergies> bestSeqsByLower;
	private final BestSet<SeqFreeEnergies> bestSeqsByUpper;

	private AffinityDirector(MultiStateConfSpace confSpace, MultiStateConfSpace.State complex, MultiStateConfSpace.State design, MultiStateConfSpace.State target, int K, double gWidthMax, Timing timing) {

		this.confSpace = confSpace;
		this.complex = complex;
		this.design = design;
		this.target = target;
		this.K = K;
		this.gWidthMax = gWidthMax;
		this.timing = timing;

		bestSeqs = new ArrayList<>();
		bestSeqsByLower = new BestSet<>(K, seq -> seq.freeEnergies[complex.index].lower);
		bestSeqsByUpper = new BestSet<>(K, seq -> seq.freeEnergies[complex.index].upper);
	}

	public void init(Directions directions, NodeProcessor processor) {

		// set the node trees to the whole space
		directions.setTrees(confSpace.states.stream()
			.map(state -> new RCs(state.confSpace))
			.toArray(RCs[]::new)
		);
	}

	public void direct(Directions directions, NodeProcessor processor) {

		directions.member.log("Searching for low free energy sequences in state %s", complex.name);
		Stopwatch stopwatch = new Stopwatch().start();

		var gcalc = new FreeEnergyCalculator();
		var ignoredSequences = new HashSet<>();
		var finishedSequences = new HashSet<Sequence>();
		var newlyFinishedSequences = new HashSet<Sequence>();

		// focus on the complex state
		directions.focus(complex.index);

		// TODO: add fail-safe to exit early if we somehow ran out of nodes to process
		while (true) {

			// wait a bit before checking progress
			ThreadTools.sleep(timing.workMs(stopwatch), TimeUnit.MILLISECONDS);

			// get all the sequences while the seqdb is locked
			Map<Sequence,StateZ[]> seqs = new HashMap<>();
			synchronized (processor.seqdb) {
				for (var entry : processor.seqdb.getSequenced()) {
					Sequence seq = entry.getKey();

					// but skip sequences we've already proven aren't in the best K
					// TODO: poplate this set?
					if (ignoredSequences.contains(seq)) {
						continue;
					}

					SeqInfo seqInfo = entry.getValue();
					seqs.put(seq, seqInfo.statezs);
				}
			}

			// compute all the free energies
			Map<Sequence,SeqFreeEnergies> seqgs = new HashMap<>();
			for (var entry : seqs.entrySet()) {
				Sequence seq = entry.getKey();
				StateZ[] statezs = entry.getValue();

				var seqg = new SeqFreeEnergies(seq, Arrays.stream(statezs)
					.map(statez -> gcalc.calc(statez.zSumBounds))
					.toArray(DoubleBounds[]::new)
				);
				seqgs.put(seq, seqg);
			}

			// mark finished sequences
			finishedSequences.clear();
			for (var seqg : seqgs.values()) {
				var complexG = seqg.freeEnergies[complex.index];
				var complexZ = seqs.get(seqg.seq)[complex.index].zSumBounds;

				// finish the sequence if the Z upper bound is negative
				// but for fully assigned sequences only
				// partially-assigned sequences often have residuals that are relatively close to zero, but sometimes slightly negative
				// (indicates not enough seqdb precision to compute the bound correctly)
				if (MathTools.isNegative(complexZ.upper) && seqg.seq.isFullyAssigned()) {
					finishedSequences.add(seqg.seq);
				}

				// finish the sequnce if the precision is on target
				if (complexG.size() <= gWidthMax) {
					finishedSequences.add(seqg.seq);
				}
			}
			directions.finishSequences(complex.sequencedIndex, finishedSequences, newlyFinishedSequences);

			// report the newly finished sequences
			for (var newlyFinishedSeq : newlyFinishedSequences) {
				var seqg = seqgs.get(newlyFinishedSeq);
				var g = seqg.freeEnergies[complex.index];
				if (Double.isNaN(g.lower)) {
					directions.member.log("WARNING: SeqDB lacks precision to compute free energy lower bound correctly for: [%s]", newlyFinishedSeq);
				}
				directions.member.log("Finished sequence: [%s]   complex G %s   width %.6f   time %s",
					newlyFinishedSeq, g, g.size(), stopwatch.getTime(2)
				);
			}

			// TODO: apply stability filters

			bestSeqs.clear();

			// find the best sequences by complex free energy bounds
			bestSeqsByLower.clear();
			bestSeqsByUpper.clear();
			for (var seqg : seqgs.values()) {
				bestSeqsByLower.add(seqg);
				bestSeqsByUpper.add(seqg);
			}

			// find all the sequences that overlap the worst best upper bound
			var worstBestUpper = bestSeqsByUpper.maxScore();
			var overlapSeqgs = new ArrayList<SeqFreeEnergies>();
			for (var seqg : seqgs.values()) {
				if (seqg.freeEnergies[complex.index].contains(worstBestUpper)) {
					overlapSeqgs.add(seqg);
				}
			}

			// are we done yet?

			// to be done, all the best sequences must have finite bounds
			if (!bestSeqsByUpper.toList().stream()
				.allMatch(seqg -> seqg.freeEnergies[complex.index].isFinite())
			) {
				reportProgress(directions, processor, finishedSequences, stopwatch);
				continue;
			}

			// to be done, all the best sequences must be fully assigned
			if (!bestSeqsByUpper.toList().stream()
				.allMatch(seqg -> seqg.seq.isFullyAssigned())
			) {
				reportProgress(directions, processor, finishedSequences, stopwatch);
				continue;
			}

			// if all of overlaps are in the best K by upper, we're done
			if (overlapSeqgs.stream()
				.allMatch(seqg -> bestSeqsByUpper.contains(seqg))
			) {
				directions.member.log("Found best %d sequences!", K);

				// which means the best K by upper and the best K by lower are the same set
				var map = new HashMap<Sequence,SeqFreeEnergies>();
				for (var seqg : bestSeqsByLower.toList()) {
					map.put(seqg.seq, seqg);
				}
				bestSeqs.addAll(map.values());

				break;
			}

			// if all the overlaps are at least the desired precision, we're done
			if (overlapSeqgs.stream()
				.allMatch(seqg -> seqg.freeEnergies[complex.index].size() <= gWidthMax)
			) {
				directions.member.log("Can't strictly find best %d sequences, but all borderline candidates are at the desired precision, so no further progress is possible.", K);

				// in this case, the best sequences are the best K by upper, and the overlaps
				var map = new HashMap<Sequence,SeqFreeEnergies>();
				for (var seqg : bestSeqsByUpper.toList()) {
					map.put(seqg.seq, seqg);
				}
				for (var seqg : overlapSeqgs) {
					map.put(seqg.seq, seqg);
				}
				bestSeqs.addAll(map.values());

				break;
			}

			// nope, not done yet. write out the most promising sequences and keep going
			reportProgress(directions, processor, finishedSequences, stopwatch);
		}

		// TODO: design
		// TODO: target

		// all done, stop the computation
		directions.stop();

		// put the best sequences in some reasonable order
		bestSeqs.sort(Comparator.comparing(seqg -> seqg.freeEnergies[complex.index].lower));

		directions.member.log(report());
	}

	private void reportProgress(Directions directions, NodeProcessor processor, Set<Sequence> finishedSeqs, Stopwatch stopwatch) {

		// show best sequences that are still unfinished
		for (var score : bestSeqsByLower.scores()) {
			var seqgs = bestSeqsByLower.values(score);

			// find the unfinished sequences
			var unfinishedSeqs = seqgs.stream()
				.filter(seqg -> !finishedSeqs.contains(seqg.seq))
				.collect(Collectors.toList());

			if (!unfinishedSeqs.isEmpty()) {

				for (var seqg : unfinishedSeqs) {
					var g = seqg.freeEnergies[complex.index];
					directions.member.log("Best unfinished sequence(s): [%s] complex G %s   width %.6f   finishedSeqs %d   nodedb %5.1f%%   time %s",
						seqg.seq,
						g.toString(3), g.size(),
						finishedSeqs.size(),
						processor.nodedb.usage()*100f,
						stopwatch.getTime(2)
					);
				}

				break;
			}
		}
	}

	public String report() {
		var buf = new StringBuilder();
		buf.append(String.format("Affinity: %d of %d sequences:\n", bestSeqsByLower.size(), K));
		for (var seqg : bestSeqs) {
			var complexG = seqg.freeEnergies[complex.index];
			var designG = seqg.freeEnergies[design.index];
			buf.append(String.format("\t[%s]   complex G %s (w %.6f)   design G %s (w %.6f)\n",
				seqg.seq,
				complexG, complexG.size(),
				designG, designG.size()
			));
		}
		// TODO: target?
		return buf.toString();
	}
}
