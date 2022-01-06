package edu.duke.cs.osprey.coffee.directors;

import edu.duke.cs.osprey.coffee.Coffee;
import edu.duke.cs.osprey.coffee.NodeProcessor;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.coffee.seqdb.SeqFreeEnergies;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;

import java.io.File;
import java.util.*;
import java.util.concurrent.TimeUnit;


public class PfuncsDirector implements Coffee.Director {

	public static class Builder {

		public final MultiStateConfSpace confSpace;
		public final MultiStateConfSpace.State state;

		private final List<Sequence> seqs = new ArrayList<>();
		private double gWidthMax = 1.0;
		private Timing timing = Timing.Efficient;
		private boolean reportProgress = false;

		//from PfuncDirector
		private int ensembleSize = 0;
		private File ensembleFile = null;
		private long ensembleUpdate = 30;
		private TimeUnit ensembleUpdateUnit = TimeUnit.SECONDS;

		public Builder(MultiStateConfSpace confSpace, MultiStateConfSpace.State state) {

			if (!state.isSequenced) {
				throw new IllegalArgumentException("Sequenced state required");
			}

			this.confSpace = confSpace;
			this.state = state;
		}

		public Builder(MultiStateConfSpace confSpace, String stateName) {
			this(confSpace, confSpace.getState(stateName));
		}

		public Builder(ConfSpace state) {
			this(
					new MultiStateConfSpace
							.Builder("state", state)
							.build(),
					"state"
			);
		}
		public Builder addSequence(Sequence val) {
			seqs.add(val);
			return this;
		}

		public Builder addSequences(Collection<Sequence> val) {
			seqs.addAll(val);
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

		public Builder setReportProgress(boolean val) {
			reportProgress = val;
			return this;
		}

		//from PfuncDirector
		/**
		 * Tracks the K lowest-energy conformations and periodically writes out an ensemble PDB file.
		 */
		public PfuncsDirector.Builder setEnsembleTracking(int size, File file) {
			ensembleSize = size;
			ensembleFile = file;
			return this;
		}

		/**
		 * Sets the minimum interval for writing the next lowest-energy ensemble PDB file.
		 */
		public PfuncsDirector.Builder setEnsembleMinUpdate(long update, TimeUnit updateUnit) {
			ensembleUpdate = update;
			ensembleUpdateUnit = updateUnit;
			return this;
		}


		public PfuncsDirector build() {
			//return new PfuncsDirector(confSpace, state, seqs, gWidthMax, timing, reportProgress, ensembleSize, ensembleFile, ensembleUpdate, ensembleUpdateUnit);
			//return new PfuncsDirector(confSpace, state, confSpace.seqSpace.getSequences(), gWidthMax, timing, reportProgress, ensembleSize, ensembleFile, ensembleUpdate, ensembleUpdateUnit);
			return new PfuncsDirector(confSpace, state, seqs, gWidthMax, timing, reportProgress, ensembleSize, ensembleFile, ensembleUpdate, ensembleUpdateUnit);



			//List<Sequence> sequencesTemp = confSpace.seqSpace.getSequences();
			//List<Sequence> inputTemp = new ArrayList<>();
			//inputTemp.add(sequencesTemp.get(0));
			//inputTemp.add(sequencesTemp.get(1));
			//return new PfuncsDirector(confSpace, state, inputTemp, gWidthMax, timing, reportProgress, ensembleSize, ensembleFile, ensembleUpdate, ensembleUpdateUnit);

		}
	}


	public final MultiStateConfSpace confSpace;
	public final MultiStateConfSpace.State state;
	public final List<Sequence> seqs;
	public final double gWidthMax;
	public final Timing timing;
	public final boolean reportProgress;
	public final List<SeqFreeEnergies> seqgs = new ArrayList<>();

	// from PfuncDirector
	public final int ensembleSize;
	public final File ensembleFile;
	public final long ensembleUpdate;
	public final TimeUnit ensembleUpdateUnit;

	private PfuncsDirector(MultiStateConfSpace confSpace, MultiStateConfSpace.State state, List<Sequence> seqs, double gWidthMax, Timing timing, boolean reportProgress, int ensembleSize, File ensembleFile, long ensembleUpdate, TimeUnit ensembleUpdateUnit) {
		this.confSpace = confSpace;
		this.state = state;
		this.seqs = seqs;
		this.gWidthMax = gWidthMax;
		this.timing = timing;
		this.reportProgress = reportProgress;
		//from PfuncDirector
		this.ensembleSize = ensembleSize;
		this.ensembleFile = ensembleFile;
		this.ensembleUpdate = ensembleUpdate;
		this.ensembleUpdateUnit = ensembleUpdateUnit;
	}

	@Override
	public void direct(Directions directions, NodeProcessor processor) {

		for (var seq : seqs) {
			seqgs.add(new SeqFreeEnergies(seq, processor.nodedb.confSpace.sequencedStates.stream()
				.map(state -> {

					if (state != this.state) {
						return null;
					}

					var pfunc = new PfuncDirector.Builder(confSpace, state)
						.setSequence(seq)
						.setGWidthMax(gWidthMax)
						.setTiming(timing)
						.setReportProgress(reportProgress)
//						.setEnsembleTracking(ensembleSize, ensembleFile)
//						.setEnsembleMinUpdate(ensembleUpdate, ensembleUpdateUnit)
						.build();
					var g = pfunc.calc(directions, processor);

					directions.finishSequence(state.sequencedIndex, seq);

					return g;
				})
				.toArray(DoubleBounds[]::new)
			));
		}

		// cleanup the nodedb when finised
		processor.nodedb.clear(state.index);
	}

}
