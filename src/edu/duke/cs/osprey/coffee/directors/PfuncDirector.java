package edu.duke.cs.osprey.coffee.directors;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.Coffee;
import edu.duke.cs.osprey.coffee.FreeEnergyCalculator;
import edu.duke.cs.osprey.coffee.NodeProcessor;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.coffee.seqdb.StateZ;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.parallelism.ThreadTools;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;


public class PfuncDirector implements Coffee.Director {

	// my kingdom for default-able arguments!

	public static class Builder {

		public final MultiStateConfSpace confSpace;
		public final MultiStateConfSpace.State state;
		public final Sequence seq;

		private double gWidthMax = 1.0;
		private Timing timing = Timing.Efficient;
		private boolean reportProgress = false;
		private int ensembleSize = 0;
		private File ensembleFile = null;
		private long ensembleUpdate = 30;
		private TimeUnit ensembleUpdateUnit = TimeUnit.SECONDS;

		public Builder(MultiStateConfSpace confSpace, MultiStateConfSpace.State state, Sequence seq) {
			this.confSpace = confSpace;
			this.state = state;
			this.seq = seq;
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

		public Builder setReportProgress(boolean val) {
			reportProgress = val;
			return this;
		}

		/**
		 * Tracks the K lowest-energy conformations and periodically writes out an ensemble PDB file.
		 */
		public Builder setEnsembleTracking(int size, File file) {
			ensembleSize = size;
			ensembleFile = file;
			return this;
		}

		/**
		 * Sets the minimum interval for writing the next lowest-energy ensemble PDB file.
		 */
		public Builder setEnsembleMinUpdate(long update, TimeUnit updateUnit) {
			ensembleUpdate = update;
			ensembleUpdateUnit = updateUnit;
			return this;
		}

		public PfuncDirector build() {
			return new PfuncDirector(confSpace, state, seq, gWidthMax, timing, reportProgress, ensembleSize, ensembleFile, ensembleUpdate, ensembleUpdateUnit);
		}
	}

	public enum Timing {

		/**
		 * Uniform 1 s batches of computation.
		 * Good for benchmarking when precise timing on short timescales is necessary.
		 */
		Precise {
			@Override
			public long workMs(Stopwatch stopwatch) {
				return 1000; // 1 s work batches
			}
		},

		/**
		 * Within the first 5 minutes of computation, use quicker batches, but gradually grow to 1 min batches.
		 * Good for long-running computations when we want to minimize overhead.
		 */
		Efficient {
			@Override
			public long workMs(Stopwatch stopwatch) {
				return Math.max(500L, (long)(Math.min(stopwatch.getTimeM()/5.0, 1.0)*60_000));
			}
		};

		public abstract long workMs(Stopwatch stopwatch);
	}

	public final MultiStateConfSpace confSpace;
	public final MultiStateConfSpace.State state;
	public final Sequence seq;
	public final double gWidthMax;
	public final Timing timing;
	public final boolean reportProgress;
	public final int ensembleSize;
	public final File ensembleFile;
	public final long ensembleUpdate;
	public final TimeUnit ensembleUpdateUnit;

	public boolean showBoundStats = false;

	private final FreeEnergyCalculator gcalc = new FreeEnergyCalculator();

	private DoubleBounds freeEnergy;

	private PfuncDirector(MultiStateConfSpace confSpace, MultiStateConfSpace.State state, Sequence seq, double gWidthMax, Timing timing, boolean reportProgress, int ensembleSize, File ensembleFile, long ensembleUpdate, TimeUnit ensembleUpdateUnit) {
		this.confSpace = confSpace;
		this.state = state;
		this.seq = seq;
		this.gWidthMax = gWidthMax;
		this.timing = timing;
		this.reportProgress = reportProgress;
		this.ensembleSize = ensembleSize;
		this.ensembleFile = ensembleFile;
		this.ensembleUpdate = ensembleUpdate;
		this.ensembleUpdateUnit = ensembleUpdateUnit;
	}

	@Override
	public int numBestConfs() {
		return ensembleSize;
	}

	@Override
	public void init(Directions directions, NodeProcessor processor) {

		// set the node trees for each state to just the specified sequence
		directions.setTrees(confSpace.states.stream()
			.map(state -> seq.makeRCs(state.confSpace))
			.toArray(RCs[]::new)
		);
	}

	@Override
	public void direct(Directions directions, NodeProcessor processor) {

		// tell the cluster to focus on this state
		directions.focus(state.index);

		freeEnergy = calc(directions, processor);

		// all done, stop the computation
		directions.stop();
	}

	DoubleBounds calc(Directions directions, NodeProcessor processor) {

		directions.member.log("Processing state %s", state.name);
		Stopwatch stopwatch = new Stopwatch().start();

		long lastEnsembleNs = stopwatch.getTimeNs();

		// wait until the free energy is precise enough
		// TODO: add fail-safe to exit early if we somehow ran out of nodes to process
		while (true) {

			// wait a bit before checking progress
			ThreadTools.sleep(timing.workMs(stopwatch), TimeUnit.MILLISECONDS);

			// lock the seqdb and get the statez
			StateZ statez;
			synchronized (processor.seqdb) {
				statez = processor.seqdb.get(state, seq);
			}

			// compute the current bounds on free energy
			DoubleBounds g = gcalc.calc(statez.zSumBounds);
			double gWidth = g.size();

			// what's the best precision we could ever get given the nodes we've already dropped?
			double gWidthMin = Math.abs(gcalc.calc(
				processor.seqdb.bigMath()
					.set(statez.zSumBounds.upper)
					.sub(statez.zSumDropped)
					.recip()
					.mult(statez.zSumBounds.upper)
					.get()
			));

			/* PROOF:
				At step i, zSum_i is bounded by [min_i, max_i].
				Let U_i = max_i - min_i, be the uncertainty at step i.
				Let D_i be the part of the uncertainty due to dropped nodes at step i.

				When the computations finishes at step n, the same still holds:
					zSum_n is bounded by [min_n, max_n].
				Assuming we didn't drop any more nodes since step i, then:
					U_n = D_i
					min_n >= min_i
					max_n <= max_i.

				By definition:
					min_n = max_n - U_n.
				Then after some substitutions:
					min_n = max_n - D_i
						 <= max_i - D_i.

				So min_n is bounded by [min_i, max_i - D_i].

				At step n, zSum_n is bounded by [min_n, min_n + U_n] by definition
					and therefore by [min_n, min_n + D_i], snce U_n = D_i.

				After converting bounds on zSum_n to free energy, the width of the free energy bounds depends
				not only on D_i, but also on the value of min_n.
				Since we're intersted in the minimim free energy width possible after dropping nodes,
				we must choose the min_n in [min_i, max_i - D_i] that minimizes the free energy width.
				The optimial value for min_n happens to be max_i - D_i, so that's convenient. =P
				(proof of the last step omitted, it's not the hard step in this proof)
			*/

			// report progress if needed
			if (reportProgress) {
				directions.member.log("\tG %s   width %.6f of %.6f   nodedb %5.1f%%   time %s",
					g.toString(3), gWidth, gWidthMin,
					processor.nodedb.usage()*100f,
					stopwatch.getTime(2)
				);
				/* TODO: show node processing statistics? eg:
					processor.stateInfos[state.index].energyBoundStats.count()
					processor.stateInfos[state.index].energyBoundStats.meanGap()
					but these are only for the local cluster member
					in general, we'd need to aggregate over all the cluster members
					maybe using some kind of statistics broadcast mechanism
				*/
				if (showBoundStats) {
					directions.member.log("%s", processor.stateInfos[state.index].energyBoundStats.toString());
				}
			}

			// should we save an ensemble now?
			if (ensembleUpdate > 0 && stopwatch.getTimeNs() >= lastEnsembleNs + ensembleUpdateUnit.toNanos(ensembleUpdate)) {
				saveEnsemble(directions, processor);
				lastEnsembleNs = stopwatch.getTimeNs();
			}

			// are we there yet?
			if (gWidth <= gWidthMax) {
				break;
			}
		}

		// save the final ensemble at the end
		saveEnsemble(directions, processor);

		return gcalc.calc(processor.seqdb.get(state, seq).zSumBounds);
	}

	public DoubleBounds getFreeEnergy() {
		return freeEnergy;
	}

	private void saveEnsemble(Directions directions, NodeProcessor processor) {

		// skip if ensemble saving is disabled
		if (ensembleFile == null) {
			return;
		}

		// get the confs
		List<int[]> bestConfs;
		synchronized (processor.seqdb) {
			bestConfs = processor.seqdb.getBestConfs(state, seq).stream()
				.map(econf -> econf.getAssignments())
				.collect(Collectors.toList());
		}

		// minimize them
		// TODO: how to handle parallelism here?
		var energiedCoords = processor.minimizeCoords(state.index, bestConfs);

		// write the PDB file
		String comment = String.format("Ensemble of %d conformations for:\n\t   State  %s\n\tSequence  [%s]",
			energiedCoords.size(), state.name, seq.toString(Sequence.Renderer.AssignmentMutations)
		);
		PDBIO.writeFileEcoords(energiedCoords, ensembleFile, comment);

		directions.member.log("Saved ensemble of %s conformations energies in [%.2f,%.2f] to %s",
			energiedCoords.size(),
			energiedCoords.stream().mapToDouble(e -> e.energy).min().orElse(Double.NaN),
			energiedCoords.stream().mapToDouble(e -> e.energy).max().orElse(Double.NaN),
			ensembleFile.getAbsolutePath()
		);
	}
}
