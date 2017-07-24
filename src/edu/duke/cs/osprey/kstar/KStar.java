package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfSearchFactory;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.math.BigDecimal;
import java.util.*;
import java.util.stream.Collectors;

public class KStar {

	public static interface ConfEnergyCalculatorFactory {
		ConfEnergyCalculator make(SimpleConfSpace confSpace, EnergyCalculator ecalc);
	}

	// *sigh* Java makes this stuff so verbose to do...
	// Kotlin would make this so much easier
	public static class Settings {

		public static class Builder {

			private double epsilon = 0.683;
			private int maxSimultaneousMutations = 1;
			private KStarScoreWriter.Writers scoreWriters = new KStarScoreWriter.Writers();

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
				return addScoreConsoleWriter(new KStarScoreWriter.Formatter.SequencePfuncsScore());
			}

			public Builder addScoreFileWriter(File file, KStarScoreWriter.Formatter val) {
				return addScoreWriter(new KStarScoreWriter.ToFile(file, val));
			}

			public Builder addScoreFileWriter(File file) {
				return addScoreFileWriter(file, new KStarScoreWriter.Formatter.Log());
			}

			public Settings build() {
				return new Settings(epsilon, maxSimultaneousMutations, scoreWriters);
			}
		}

		public final double epsilon;
		public final int maxSimultaneousMutations;
		public final KStarScoreWriter.Writers scoreWriters;

		public Settings(double epsilon, int maxSimultaneousMutations, KStarScoreWriter.Writers scoreWriters) {
			this.epsilon = epsilon;
			this.maxSimultaneousMutations = maxSimultaneousMutations;
			this.scoreWriters = scoreWriters;
		}
	}

	public static class Sequence extends ArrayList<String> {

		public Sequence(int size) {
			super(size);
		}

		public Sequence(String resTypes) {
			this(Arrays.asList(resTypes.split(" ")));
		}

		public Sequence(List<String> resTypes) {
			super(resTypes);
		}

		@Override
		public String toString() {
			return String.join(" ", this);
		}

		public String toString(List<SimpleConfSpace.Position> positions) {
			return String.join(" ", positions.stream()
				.map((pos) -> this.get(pos.index) + "-" + pos.resNum)
				.collect(Collectors.toList())
			);
		}
	}

	public static enum ConfSpaceType {
		Protein,
		Ligand,
		Complex
	}

	public class ConfSpaceInfo {

		public final ConfSpaceType type;
		public final SimpleConfSpace confSpace;
		public final ConfEnergyCalculator confEcalc;

		public final List<Sequence> sequences = new ArrayList<>();
		public EnergyMatrix emat = null;
		public final Map<Sequence,PartitionFunction.Result> pfuncResults = new HashMap<>();

		public ConfSpaceInfo(ConfSpaceType type, SimpleConfSpace confSpace) {

			this.type = type;
			this.confSpace = confSpace;
			this.confEcalc = confEcalcFactory.make(confSpace, ecalc);
		}

		public void calcEmat() {
			emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.build()
				.calcEnergyMatrix();
		}

		public Sequence makeWildTypeSequence() {
			Sequence sequence = new Sequence(confSpace.positions.size());
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				sequence.add(pos.strand.mol.getResByPDBResNumber(pos.resNum).template.name);
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

			// prune down to just this sequence
			// TODO: inherit any pruning from e.g. DEE
			PruningMatrix pmat = new PruningMatrix(confSpace, 0.0);
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				String resType = sequence.get(pos.index);

				for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {

					// if this RC doesn't match the sequence, prune it
					if (!rc.template.name.equals(resType)) {
						pmat.setOneBody(pos.index, rc.index, true);
					}
				}
			}

			// make the partition function
			PartitionFunction pfunc = new SimplePartitionFunction(emat, pmat, confSearchFactory, confEcalc);
			pfunc.setReportProgress(false); // TODO: expose option to turn off

			// compute it
			pfunc.init(settings.epsilon);
			pfunc.compute();

			// save the result
			result = pfunc.makeResult();
			pfuncResults.put(sequence, result);
			return result;
		}
	}

	public final Strand protein;
	public final Strand ligand;
	public final EnergyCalculator ecalc;
	public final ConfEnergyCalculatorFactory confEcalcFactory;
	public final ConfSearchFactory confSearchFactory;
	public final Settings settings;

	public final ConfSpaceInfo proteinInfo;
	public final ConfSpaceInfo ligandInfo;
	public final ConfSpaceInfo complexInfo;

	public final List<BigDecimal> kstarScores;

	public KStar(Strand protein, Strand ligand, SimpleConfSpace complexConfSpace, EnergyCalculator ecalc, ConfEnergyCalculatorFactory confEcalcFactory, ConfSearchFactory confSearchFactory, Settings settings) {

		this.protein = protein;
		this.ligand = ligand;
		this.ecalc = ecalc;
		this.confEcalcFactory = confEcalcFactory;
		this.confSearchFactory = confSearchFactory;
		this.settings = settings;

		// make conf space infos
		proteinInfo = new ConfSpaceInfo(ConfSpaceType.Protein, complexConfSpace.makeSubspace(protein));
		ligandInfo = new ConfSpaceInfo(ConfSpaceType.Ligand, complexConfSpace.makeSubspace(ligand));
		complexInfo = new ConfSpaceInfo(ConfSpaceType.Complex, complexConfSpace);

		kstarScores = new ArrayList<>();
	}

	public void run() {

		// compute energy matrices
		proteinInfo.calcEmat();
		ligandInfo.calcEmat();
		complexInfo.calcEmat();

		Map<SimpleConfSpace.Position,SimpleConfSpace.Position> complexToProteinMap = complexInfo.confSpace.mapPositionsTo(proteinInfo.confSpace);
		Map<SimpleConfSpace.Position,SimpleConfSpace.Position> complexToLigandMap = complexInfo.confSpace.mapPositionsTo(ligandInfo.confSpace);

		// collect the wild type sequences
		proteinInfo.sequences.add(proteinInfo.makeWildTypeSequence());
		ligandInfo.sequences.add(ligandInfo.makeWildTypeSequence());
		complexInfo.sequences.add(complexInfo.makeWildTypeSequence());

		// collect all the sequences explicitly
		List<List<SimpleConfSpace.Position>> powersetOfPositions = MathTools.powersetUpTo(complexInfo.confSpace.positions, settings.maxSimultaneousMutations);
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
				Sequence complexSequence = complexInfo.makeWildTypeSequence();
				for (int i=0; i<mutablePositions.size(); i++) {
					complexSequence.set(mutablePositions.get(i).index, mutations.get(i));
				}
				complexInfo.sequences.add(complexSequence);

				// split complex sequence into protein/ligand sequences
				Sequence proteinSequence = proteinInfo.makeWildTypeSequence();
				Sequence ligandSequence = ligandInfo.makeWildTypeSequence();
				for (SimpleConfSpace.Position pos : complexInfo.confSpace.positions) {

					SimpleConfSpace.Position proteinPos = complexToProteinMap.get(pos);
					if (proteinPos != null) {
						proteinSequence.set(proteinPos.index, complexSequence.get(pos.index));
					}

					SimpleConfSpace.Position ligandPos = complexToLigandMap.get(pos);
					if (ligandPos != null) {
						ligandSequence.set(ligandPos.index, complexSequence.get(pos.index));
					}
				}
				proteinInfo.sequences.add(proteinSequence);
				ligandInfo.sequences.add(ligandSequence);
			}
		}

		// TODO: sequence filtering? do we need to reject some mutation combinations for some reason?

		System.out.println("computing K* scores for " + complexInfo.sequences.size() + " sequences...");
		settings.scoreWriters.writeHeader();
		// TODO: progress bar?

		// compute all the partition functions and K* scores
		int n = complexInfo.sequences.size();
		for (int i=0; i<n; i++) {

			// get the pfuncs
			PartitionFunction.Result proteinResult = proteinInfo.calcPfunc(i);
			PartitionFunction.Result ligandResult = ligandInfo.calcPfunc(i);
			PartitionFunction.Result complexResult = complexInfo.calcPfunc(i);

			// compute the K* score if possible
			BigDecimal kstarScore = PartitionFunction.Result.calcKStarScore(proteinResult, ligandResult, complexResult);
			kstarScores.add(kstarScore);

			// report scores
			settings.scoreWriters.writeScore(new KStarScoreWriter.ScoreInfo(
				i,
				n,
				complexInfo.sequences.get(i),
				complexInfo.confSpace,
				proteinResult, ligandResult, complexResult,
				kstarScore
			));
		}
	}
}
