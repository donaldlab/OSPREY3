package edu.duke.cs.osprey.coffee.directors;

import edu.duke.cs.osprey.coffee.Coffee;
import edu.duke.cs.osprey.coffee.NodeProcessor;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.coffee.seqdb.SeqFreeEnergies;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;


/**
 * An implementation of the K* algorithm that uses
 * COFFEE for partition function calculation.
 *
 * The complex and design should be mutable, the target should not.
 */
public class KStarDirector implements Coffee.Director {

	public static class Builder {

		private final MultiStateConfSpace confSpace;
		private final MultiStateConfSpace.State complex;
		private final MultiStateConfSpace.State design;
		private final MultiStateConfSpace.State target;

		/**
		 * Sets the largest precision desired for free energy calculations.
		 * Resulting free energy values may be more precise than this maximum value.
		 * If the computation is stopped early, the free energy values may be less precise than this value.
		 */
		private double gWidthMax = 1.0;

		/**
		 * Pruning criteria to remove sequences with unstable unbound states relative to the wild type sequence.
		 * Defined in units of kcal/mol.
		 *
		 * Set to null to disable the filter entirely.
		 *
		 * If the the wild-type sequence is not in the list of sequences to compute,
		 * and the stability threshold is not null,
		 * then the free energies for the wild-type sequence will be calculated anyway.
		 *
		 * More precisely, a sequence is pruned when the following expression is true:
		 *
		 * L(G_s) > U(G_w) + t
		 *
		 * where:
		 *   - L(G_s) is the lower bound on the free energy for sequence s
		 *   - U(G_w) is the upper bound on the free energy for the wild-type sequence w
		 *   - t is the stability threshold
		 */
		private Double stabilityThreshold = 5.0;

		private Integer maxSimultaneousMutations = 1;

		private Timing timing = Timing.Efficient;

		/**
		 * Set to true to report progress on computing individual conf space states to the log.
		 */
		private boolean reportStateProgress = true;

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

		public Builder setGWidthMax(double val) {
			gWidthMax = val;
			return this;
		}

		public Builder setMaxSimultaneousMutations(Integer val) {
			maxSimultaneousMutations = val;
			return this;
		}

		public Builder setStabilityThreshold(Double val) {
			stabilityThreshold = val;
			return this;
		}

		public Builder setTiming(Timing val) {
			timing = val;
			return this;
		}

		public Builder setReportStateProgress(boolean val) {
			reportStateProgress = val;
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

		public KStarDirector build() {

			// check for common errors early and give them friendlier error messages
			if (stabilityThreshold != null && !confSpace.seqSpace.containsWildTypeSequence()) {
				throw new IllegalArgumentException(
					"Stability threshold is enabled, but the sequence space does not define the wild-type sequence."
					+ "\nEither disable the stability threshold, or use a sequence space that contains the wild-type sequence."
				);
			}

			return new KStarDirector(
				confSpace, complex, design, target,
				gWidthMax, stabilityThreshold, maxSimultaneousMutations, timing, reportStateProgress,
				ensembleSize, ensembleDir, ensembleUpdate, ensembleUpdateUnit
			);
		}
	}

	public final MultiStateConfSpace confSpace;
	public final MultiStateConfSpace.State complex;
	public final MultiStateConfSpace.State design;
	public final MultiStateConfSpace.State target;
	public final double gWidthMax;
	public final Double stabilityThreshold;
	public final Integer maxSimultaneousMutations;
	public final Timing timing;
	public final boolean reportStateProgress;
	public final int ensembleSize;
	public final File ensembleDir;
	public final long ensembleUpdate;
	public final TimeUnit ensembleUpdateUnit;

	public SeqFreeEnergies wildTypeG = null;
	public final List<Sequence> sequences = new ArrayList<>();
	public final List<SeqFreeEnergies> sequenceGs = new ArrayList<>();
	public DoubleBounds targetG = null;


	private KStarDirector(
		MultiStateConfSpace confSpace, MultiStateConfSpace.State complex, MultiStateConfSpace.State design, MultiStateConfSpace.State target,
		double gWidthMax, Double stabilityThreshold, Integer maxSimultaneousMutations, Timing timing, boolean reportStateProgress,
		int ensembleSize, File ensembleDir, long ensembleUpdate, TimeUnit ensembleUpdateUnit
	) {

		this.confSpace = confSpace;
		this.complex = complex;
		this.design = design;
		this.target = target;
		this.gWidthMax = gWidthMax;
		this.stabilityThreshold = stabilityThreshold;
		this.maxSimultaneousMutations = maxSimultaneousMutations;
		this.timing = timing;
		this.reportStateProgress = reportStateProgress;
		this.ensembleSize = ensembleSize;
		this.ensembleDir = ensembleDir;
		this.ensembleUpdate = ensembleUpdate;
		this.ensembleUpdateUnit = ensembleUpdateUnit;

		if (maxSimultaneousMutations != null) {
			// collect all the sequences explicitly, because that's how K* rolls
			var seqSpace = complex.confSpace.seqSpace();
			if (seqSpace.containsWildTypeSequence()) {
				sequences.add(seqSpace.makeWildTypeSequence());
			}
			sequences.addAll(seqSpace.getMutants(maxSimultaneousMutations));
		}
	}

	@Override
	public int numBestConfs() {
		return ensembleSize;
	}

	@Override
	public void direct(Directions directions, NodeProcessor processor) {

		var stopwatch = new Stopwatch().start();
		var seqSpace = complex.confSpace.seqSpace();

		// fill the output with nulls
		sequenceGs.clear();
		for (int i=0; i<sequences.size(); i++) {
			sequenceGs.add(null);
		}

		directions.member.log("Computing free energies for %d sequences...", sequences.size());

		// did the user request the wild-type sequence?
		int wtIndex = IntStream.range(0, sequences.size())
			.filter(i -> sequences.get(i).isWildType())
			.findFirst()
			.orElse(-1);

		// compute the wild-type sequence first, if necessary
		Double gMaxDesign = null;
		if (stabilityThreshold != null || wtIndex > 0) {

			var seq = seqSpace.makeWildTypeSequence();

			var g = new DoubleBounds[confSpace.states.size()];
			wildTypeG = new SeqFreeEnergies(seq, g);
			if (wtIndex > 0) {
				sequenceGs.set(wtIndex, wildTypeG);
			}

			// compute the free energies
			directions.member.log("wild-type sequence: %s ...", sequences.size(), seq);
			g[target.index] = calcG(directions, processor, target, null, null);
			directions.member.log("\t%20s: %s", target.name, gToString(g[target.index]));
			g[design.index] = calcG(directions, processor, design, seq, null);
			directions.member.log("\t%20s: %s", design.name, gToString(g[design.index]));
			g[complex.index] = calcG(directions, processor, complex, seq, null);
			directions.member.log("\t%20s: %s", complex.name, gToString(g[complex.index]));

			// specialize the stability threshold to the design state
			if (stabilityThreshold != null) {
				gMaxDesign = g[design.index].upper + stabilityThreshold;
			}

			// save the target G for later, since it's the same for all sequences
			targetG = g[target.index];
		}

		// then compute all the other sequences serially
		for (int i=0; i<sequences.size(); i++) {

			// skip the wild-type if we already computed it
			if (i == wtIndex) {
				continue;
			}

			var seq = sequences.get(i);
			directions.member.log("%4d/%d: sequence: %s ...", i + 1, sequences.size(), seq);

			var g = new DoubleBounds[confSpace.states.size()];
			sequenceGs.set(i, new SeqFreeEnergies(seq, g));

			// compute target free energy, if we don't have it already
			if (targetG == null) {
				g[target.index] = calcG(directions, processor, target, null, null);
				directions.member.log("\t%20s: %s", target.name, gToString(g[target.index]));
				targetG = g[target.index];
			}

			// compute the design free energy
			g[design.index] = calcG(directions, processor, design, seq, gMaxDesign);
			directions.member.log("\t%20s: %s", design.name, gToString(g[design.index]));

			// check the stabiliby threshold, if needed, ie: L(G_s) > U(G_w) + t
			if (gMaxDesign != null && g[design.index].lower > gMaxDesign) {
				directions.member.log("\t%20s  Sequence pruned, failed stability filter", "");
				continue;
			}

			// finally, compute the complex free energy
			g[complex.index] = calcG(directions, processor, complex, seq, null);
			directions.member.log("\t%20s: %s", complex.name, gToString(g[complex.index]));
		}

		directions.member.log("K* calculation complete in %s!", stopwatch.getTime(2));
	}

	private DoubleBounds calcG(Directions directions, NodeProcessor processor, MultiStateConfSpace.State state, Sequence seq, Double gMax) {

		processor.nodedb.clear(state.index);

		var pfunc = new PfuncDirector.Builder(confSpace, state, seq)
			.setGWidthMax(gWidthMax)
			.setTiming(timing)
			.setReportProgress(reportStateProgress)
			.setGMax(gMax);

		// for the complex state, pass along the ensemble reporting settings
		if (state == complex) {
			pfunc.setEnsembleTracking(ensembleSize, new File(ensembleDir, String.format("seq.%s.pdb", seq.toString(Sequence.Renderer.AssignmentMutations))));
			pfunc.setEnsembleMinUpdate(ensembleUpdate, ensembleUpdateUnit);
		}

		return pfunc.build().calc(directions, processor);
	}

	private String gToString(DoubleBounds g) {
		return String.format("%s (w %.6f) kcal/mol", g.toString(6), g.size());
	}
}
