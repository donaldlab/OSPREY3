package edu.duke.cs.osprey.coffee.directors;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.Coffee;
import edu.duke.cs.osprey.coffee.FreeEnergyCalculator;
import edu.duke.cs.osprey.coffee.NodeProcessor;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.coffee.nodedb.NodeTree;
import edu.duke.cs.osprey.coffee.seqdb.SeqFreeEnergies;
import edu.duke.cs.osprey.coffee.seqdb.SeqInfo;
import edu.duke.cs.osprey.coffee.seqdb.StateZ;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.parallelism.ThreadTools;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
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
		private Integer maxSimultaneousMutations = 1;
		private Timing timing = Timing.Efficient;
		private int ensembleSize = 0;
		private File ensembleDir = null;
		private long ensembleUpdate = 5;
		private TimeUnit ensembleUpdateUnit = TimeUnit.MINUTES;


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

		public Builder(ConfSpace complex, ConfSpace design, ConfSpace target) {
			this(
					new MultiStateConfSpace
							.Builder("complex", complex)
							.addMutableState("design", design)
							.addUnmutableState("target", target)
							.build(),
					"complex",
					"design",
					"target"
			);
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

		public Builder setMaxSimultaneousMutations(Integer val) {
			maxSimultaneousMutations = val;
			return this;
		}

		public Builder setTiming(Timing val) {
			timing = val;
			return this;
		}

		/**
		 * Tracks the lowest-energy conformations for the best sequences
		 * and periodically writes out ensemble PDB files to the specified directory.
		 */
		public Builder setEnsembleTracking(int size, File dir) {
			ensembleSize = size;
			ensembleDir = dir;
			return this;
		}

		/**
		 * Sets the minimum interval for writing the next lowest-energy ensemble PDB files.
		 */
		public Builder setEnsembleMinUpdate(long update, TimeUnit updateUnit) {
			ensembleUpdate = update;
			ensembleUpdateUnit = updateUnit;
			return this;
		}

		public AffinityDirector build() {
			return new AffinityDirector(
				confSpace, complex, design, target,
				K, gWidthMax, maxSimultaneousMutations, timing,
				ensembleSize, ensembleDir, ensembleUpdate, ensembleUpdateUnit
			);
		}
	}

	public final MultiStateConfSpace confSpace;
	public final MultiStateConfSpace.State complex;
	public final MultiStateConfSpace.State design;
	public final MultiStateConfSpace.State target;
	public final int K;
	public final double gWidthMax;
	public final Integer maxSimultaneousMutations;
	public final Timing timing;
	public final int ensembleSize;
	public final File ensembleDir;
	public final long ensembleUpdate;
	public final TimeUnit ensembleUpdateUnit;

	public final List<SeqFreeEnergies> bestSeqs;
	public DoubleBounds targetFreeEnergy;

	private final BestSet<SeqFreeEnergies> bestSeqsByLower;
	private final BestSet<SeqFreeEnergies> bestSeqsByUpper;

	private AffinityDirector(MultiStateConfSpace confSpace, MultiStateConfSpace.State complex, MultiStateConfSpace.State design, MultiStateConfSpace.State target, int K, double gWidthMax, Integer maxSimultaneousMutations, Timing timing, int ensembleSize, File ensembleDir, long ensembleUpdate, TimeUnit ensembleUpdateUnit) {

		this.confSpace = confSpace;
		this.complex = complex;
		this.design = design;
		this.target = target;
		this.K = K;
		this.gWidthMax = gWidthMax;
		this.maxSimultaneousMutations = maxSimultaneousMutations;
		this.timing = timing;
		this.ensembleSize = ensembleSize;
		this.ensembleDir = ensembleDir;
		this.ensembleUpdate = ensembleUpdate;
		this.ensembleUpdateUnit = ensembleUpdateUnit;

		bestSeqs = new ArrayList<>();
		targetFreeEnergy = null;
		bestSeqsByLower = new BestSet<>(K, seq -> seq.freeEnergies[complex.index].lower);
		bestSeqsByUpper = new BestSet<>(K, seq -> seq.freeEnergies[complex.index].upper);
	}

	@Override
	public int numBestConfs() {
		return ensembleSize;
	}

	@Override
	public void direct(Directions directions, NodeProcessor processor) {

		directions.member.log("Searching for low free energy sequences in state %s", complex.name);
		Stopwatch stopwatch = new Stopwatch().start();

		bestSeqs.clear();
		targetFreeEnergy = null;

		// init the cluster with the complex state, all sequences
		directions.focus(complex.index);
		var tree = new NodeTree(new RCs(complex.confSpace), maxSimultaneousMutations);
		directions.setTree(complex.index, tree);
		processor.initRootNode(complex.index, tree);

		var gcalc = new FreeEnergyCalculator();
		var finishedSequences = new HashSet<Sequence>();
		var newlyFinishedSequences = new HashSet<Sequence>();

		long lastEnsembleNs = stopwatch.getTimeNs();

		// first, find the sequences with the best complex free energies
		directions.focus(complex.index);
		// TODO: add fail-safe to exit early if we somehow ran out of nodes to process?
		while (true) {

			// wait a bit before checking progress
			ThreadTools.sleep(timing.workMs(stopwatch), TimeUnit.MILLISECONDS);

			// get all the sequences while the seqdb is locked
			Map<Sequence,StateZ[]> seqs = new HashMap<>();
			synchronized (processor.seqdb) {
				for (var entry : processor.seqdb.getSequenced()) {
					Sequence seq = entry.getKey();

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

			// should we save ensembles now?
			if (ensembleUpdate > 0 && stopwatch.getTimeNs() >= lastEnsembleNs + ensembleUpdateUnit.toNanos(ensembleUpdate)) {
				saveEnsemble(directions, processor, bestSeqsByLower.toList());
				lastEnsembleNs = stopwatch.getTimeNs();
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
		processor.nodedb.clear(complex.index);

		// put the best sequences in some reasonable order
		bestSeqs.sort(Comparator.comparing(seqg -> seqg.freeEnergies[complex.index].lower));

		// save the last ensembles
		saveEnsemble(directions, processor, bestSeqs);

		// next, compute the design state free energies for the best sequences
		directions.member.log("Computing design state free energies for %d sequences ...", bestSeqs.size());
		for (var seqg : bestSeqs) {
			processor.nodedb.clear(design.index);
			var pfunc = new PfuncDirector.Builder(confSpace, design, seqg.seq)
				.setGWidthMax(gWidthMax)
				.setTiming(timing)
				.setReportProgress(true) // TODO: expose setting
				.build();
			seqg.freeEnergies[design.index] = pfunc.calc(directions, processor);
		}

		// finally, compute the target state free energy
		directions.member.log("Computing target state free energies ...");
		processor.nodedb.clear(target.index);
		var pfunc = new PfuncDirector.Builder(confSpace, target)
			.setGWidthMax(gWidthMax)
			.setTiming(timing)
			.setReportProgress(true) // TODO: expose setting
			.build();
		targetFreeEnergy = pfunc.calc(directions, processor);

		// all done
		directions.member.log("All three states complete in %s", stopwatch.getTime(2));
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

	private void saveEnsemble(Directions directions, NodeProcessor processor, List<SeqFreeEnergies> seqgs) {

		// skip if ensemble saving is disabled
		if (ensembleDir == null) {
			return;
		}

		// keep only fully-assigned sequences
		seqgs = seqgs.stream()
			.filter(seqg -> seqg.seq.isFullyAssigned())
			.collect(Collectors.toList());

		var stopwatch = new Stopwatch().start();
		directions.member.log("Writing ensembles for %d sequences ...", seqgs.size());
		ensembleDir.mkdirs();

		for (var seqg : seqgs) {

			// pick a name for the ensemble file
			String seqstr = seqg.seq.toString(Sequence.Renderer.ResTypeMutations)
				.replace(' ', '-');
			File ensembleFile = new File(ensembleDir, String.format("seq.%s.pdb", seqstr));

			// get the confs, if any
			List<int[]> bestConfs;
			synchronized (processor.seqdb) {
				bestConfs = processor.seqdb.getBestConfs(complex, seqg.seq).stream()
					.map(econf -> econf.getAssignments())
					.collect(Collectors.toList());
			}
			if (bestConfs.isEmpty()) {
				// no structures, write an empty file
				PDBIO.writeFileEcoords(Collections.emptyList(), ensembleFile, "No structures for this sequence");
				continue;
			}

			// minimize them
			// TODO: optimize somehow? cache structures?
			var energiedCoords = processor.minimizeCoords(complex.index, bestConfs);

			// write the PDB file
			String comment = String.format("Ensemble of %d conformations for:\n\t   State  %s\n\tSequence  [%s]",
				energiedCoords.size(), complex.name, seqg.seq.toString(Sequence.Renderer.AssignmentMutations)
			);
			PDBIO.writeFileEcoords(energiedCoords, ensembleFile, comment);
		}

		directions.member.log("Saved ensembles to %s in %s", ensembleDir.getAbsolutePath(), stopwatch.getTime(2));
	}
}
