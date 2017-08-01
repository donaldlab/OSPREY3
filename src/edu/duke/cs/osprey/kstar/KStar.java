package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public class KStar {

	public static interface ConfEnergyCalculatorFactory {
		ConfEnergyCalculator make(SimpleConfSpace confSpace, EnergyCalculator ecalc);
	}

	public static interface ConfSearchFactory {
		public ConfSearch make(EnergyMatrix emat, RCs rcs);
	}

	// *sigh* Java makes this stuff so verbose to do...
	// Kotlin would make this so much easier
	public static class Settings {

		public static class Builder {

			private double epsilon = 0.683;
			private int maxSimultaneousMutations = 1;
			private KStarScoreWriter.Writers scoreWriters = new KStarScoreWriter.Writers();
			private boolean showPfuncProgress = false;

			public Builder setEpsilon(double val) {
				epsilon = val;
				return this;
			}

			public Builder setMaxSimultaneousMutations(int val) {
				maxSimultaneousMutations = val;
				return this;
			}

			public Builder addScoreWriter(KStarScoreWriter val) {
				scoreWriters.add(val);
				return this;
			}

			public Builder addScoreConsoleWriter(KStarScoreWriter.Formatter val) {
				return addScoreWriter(new KStarScoreWriter.ToConsole(val));
			}

			public Builder addScoreConsoleWriter() {
				return addScoreConsoleWriter(new KStarScoreWriter.Formatter.SequenceKStarPfuncs());
			}

			public Builder addScoreFileWriter(File file, KStarScoreWriter.Formatter val) {
				return addScoreWriter(new KStarScoreWriter.ToFile(file, val));
			}

			public Builder addScoreFileWriter(File file) {
				return addScoreFileWriter(file, new KStarScoreWriter.Formatter.Log());
			}

			public Builder setShowPfuncProgress(boolean val) {
				showPfuncProgress = val;
				return this;
			}

			public Settings build() {
				return new Settings(epsilon, maxSimultaneousMutations, scoreWriters, showPfuncProgress);
			}
		}

		public final double epsilon;
		public final int maxSimultaneousMutations;
		public final KStarScoreWriter.Writers scoreWriters;
		public final boolean showPfuncProgress;

		public Settings(double epsilon, int maxSimultaneousMutations, KStarScoreWriter.Writers scoreWriters, boolean dumpPfuncConfs) {
			this.epsilon = epsilon;
			this.maxSimultaneousMutations = maxSimultaneousMutations;
			this.scoreWriters = scoreWriters;
			this.showPfuncProgress = dumpPfuncConfs;
		}
	}

	public static class Sequence extends ArrayList<String> {

		public static Sequence makeWildType(SimpleConfSpace confSpace) {
			return new KStar.Sequence(confSpace.positions.size()).fillWildType(confSpace);
		}

		public Sequence(int size) {
			super(size);
			for (int i=0; i<size; i++) {
				add(null);
			}
		}

		public Sequence(String resTypes) {
			this(Arrays.asList(resTypes.split(" ")));
		}

		public Sequence(List<String> resTypes) {
			super(resTypes);
		}

		public Sequence makeWithAssignment(int posIndex, String resType) {
			Sequence assigned = new Sequence(this);
			assigned.set(posIndex, resType);
			return assigned;
		}

		public Sequence fillWildType(SimpleConfSpace confSpace) {
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				if (get(pos.index) == null) {
					set(pos.index, pos.strand.mol.getResByPDBResNumber(pos.resNum).template.name);
				}
			}
			return this;
		}

		public int countAssignments() {
			int count = 0;
			for (String resType : this) {
				if (resType != null) {
					count++;
				}
			}
			return count;
		}

		public boolean isFullyAssigned() {
			for (String resType : this) {
				if (resType == null) {
					return false;
				}
			}
			return true;
		}

		public int countMutations(SimpleConfSpace confSpace) {
			int count = 0;
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				String resType = get(pos.index);
				if (resType != null && !resType.equals(pos.resFlex.wildType)) {
					count++;
				}
			}
			return count;
		}

		@Override
		public String toString() {
			return String.join(" ", this);
		}

		/** highlight mutations in upper case, show wild type residues in lower case */
		public String toString(KStar.Sequence wildtype) {
			List<String> resTypes = new ArrayList<>();
			for (int i=0; i<size(); i++) {
				String resType = get(i);
				String resWildType = wildtype.get(i);
				if (resType == null) {
					resTypes.add(null);
				} else if (resType.equals(resWildType)) {
					resTypes.add(resType.toLowerCase());
				} else {
					resTypes.add(resType.toUpperCase());
				}
			}
			return String.join(" ", resTypes);
		}

		public String toString(List<SimpleConfSpace.Position> positions) {
			return String.join(" ", positions.stream()
				.map((pos) -> this.get(pos.index) + "-" + pos.resNum)
				.collect(Collectors.toList())
			);
		}
	}

	public class ConfSpaceInfo {

		public final SimpleConfSpace confSpace;
		public final ConfEnergyCalculator confEcalc;

		public final List<Sequence> sequences = new ArrayList<>();
		public EnergyMatrix emat = null;
		public final Map<Sequence,PartitionFunction.Result> pfuncResults = new HashMap<>();

		public ConfSpaceInfo(SimpleConfSpace confSpace, ConfEnergyCalculator confEcalc) {
			this.confSpace = confSpace;
			this.confEcalc = confEcalc;
		}

		public void calcEmat() {
			emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.build()
				.calcEnergyMatrix();
		}

		public Sequence makeWildTypeSequence() {
			Sequence sequence = new Sequence(confSpace.positions.size());
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				sequence.set(pos.index, pos.strand.mol.getResByPDBResNumber(pos.resNum).template.name);
			}
			return sequence;
		}

		public PartitionFunction.Result calcPfunc(int sequenceIndex) {

			Sequence sequence = sequences.get(sequenceIndex);

			// check the cache first
			PartitionFunction.Result result = pfuncResults.get(sequence);
			if (result != null) {
				return result;
			}

			// cache miss, need to compute the partition function

			// get RCs for just this sequence
			RCs rcs = new RCs(confSpace, (pos, resConf) -> {
				return resConf.template.name.equals(sequence.get(pos.index));
			});

			// make the partition function
			ConfSearch astar = confSearchFactory.make(emat, rcs);
			PartitionFunction pfunc = new SimplePartitionFunction(astar, confEcalc);
			pfunc.setReportProgress(settings.showPfuncProgress);

			// compute it
			pfunc.init(settings.epsilon);
			pfunc.compute();

			// save the result
			result = pfunc.makeResult();
			pfuncResults.put(sequence, result);
			return result;
		}
	}

	public final ConfSpaceInfo protein;
	public final ConfSpaceInfo ligand;
	public final ConfSpaceInfo complex;
	public final ConfEnergyCalculatorFactory confEcalcFactory;
	public final ConfSearchFactory confSearchFactory;
	public final Settings settings;

	public KStar(SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex, EnergyCalculator ecalc, ConfEnergyCalculatorFactory confEcalcFactory, ConfSearchFactory confSearchFactory, Settings settings) {
		this.protein = new ConfSpaceInfo(protein, confEcalcFactory.make(protein, ecalc));
		this.ligand = new ConfSpaceInfo(ligand, confEcalcFactory.make(ligand, ecalc));
		this.complex = new ConfSpaceInfo(complex, confEcalcFactory.make(complex, ecalc));
		this.confEcalcFactory = confEcalcFactory;
		this.confSearchFactory = confSearchFactory;
		this.settings = settings;
	}

	public List<KStarScore> run() {

		List<KStarScore> scores = new ArrayList<>();

		// compute energy matrices
		protein.calcEmat();
		ligand.calcEmat();
		complex.calcEmat();

		Map<SimpleConfSpace.Position,SimpleConfSpace.Position> complexToProteinMap = complex.confSpace.mapPositionsTo(protein.confSpace);
		Map<SimpleConfSpace.Position,SimpleConfSpace.Position> complexToLigandMap = complex.confSpace.mapPositionsTo(ligand.confSpace);

		// collect the wild type sequences
		protein.sequences.add(protein.makeWildTypeSequence());
		ligand.sequences.add(ligand.makeWildTypeSequence());
		complex.sequences.add(complex.makeWildTypeSequence());

		// collect all the sequences explicitly
		List<List<SimpleConfSpace.Position>> powersetOfPositions = MathTools.powersetUpTo(complex.confSpace.positions, settings.maxSimultaneousMutations);
		Collections.reverse(powersetOfPositions); // NOTE: reverse to match order of old code
		for (List<SimpleConfSpace.Position> mutablePositions : powersetOfPositions) {

			// collect the mutations (res types except for wild type) for these positions into a simple list list
			List<List<String>> resTypes = new ArrayList<>();
			for (SimpleConfSpace.Position pos : mutablePositions) {
				resTypes.add(pos.resFlex.resTypes.stream()
					.filter((resType) -> !resType.equals(pos.resFlex.wildType))
					.collect(Collectors.toList())
				);
			}

			// enumerate all the combinations of res types
			for (List<String> mutations : MathTools.cartesianProduct(resTypes)) {

				// build the complex sequence
				Sequence complexSequence = complex.makeWildTypeSequence();
				for (int i=0; i<mutablePositions.size(); i++) {
					complexSequence.set(mutablePositions.get(i).index, mutations.get(i));
				}
				complex.sequences.add(complexSequence);

				// split complex sequence into protein/ligand sequences
				Sequence proteinSequence = protein.makeWildTypeSequence();
				Sequence ligandSequence = ligand.makeWildTypeSequence();
				for (SimpleConfSpace.Position pos : complex.confSpace.positions) {

					SimpleConfSpace.Position proteinPos = complexToProteinMap.get(pos);
					if (proteinPos != null) {
						proteinSequence.set(proteinPos.index, complexSequence.get(pos.index));
					}

					SimpleConfSpace.Position ligandPos = complexToLigandMap.get(pos);
					if (ligandPos != null) {
						ligandSequence.set(ligandPos.index, complexSequence.get(pos.index));
					}
				}
				protein.sequences.add(proteinSequence);
				ligand.sequences.add(ligandSequence);
			}
		}

		// TODO: sequence filtering? do we need to reject some mutation combinations for some reason?

		System.out.println("computing K* scores for " + complex.sequences.size() + " sequences to epsilon = " + settings.epsilon + " ...");
		settings.scoreWriters.writeHeader();
		// TODO: progress bar?

		// compute all the partition functions and K* scores
		int n = complex.sequences.size();
		for (int i=0; i<n; i++) {

			// get the pfuncs
			PartitionFunction.Result proteinResult = protein.calcPfunc(i);
			PartitionFunction.Result ligandResult = ligand.calcPfunc(i);
			PartitionFunction.Result complexResult = complex.calcPfunc(i);

			// compute the K* score if possible
			KStarScore kstarScore = new KStarScore(proteinResult, ligandResult, complexResult);
			scores.add(kstarScore);

			// report scores
			settings.scoreWriters.writeScore(new KStarScoreWriter.ScoreInfo(
				i,
				n,
				complex.sequences.get(i),
				complex.confSpace,
				kstarScore
			));
		}

		return scores;
	}
}
