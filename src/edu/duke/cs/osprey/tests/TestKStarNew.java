package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStar.ConfSearchFactory;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.util.HashSet;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Function;


public class TestKStarNew {

	public static void main(String[] args) {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readFile("examples/python.KStar/2RL0.min.reduce.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.addMoleculeForWildTypeRotamers(mol)
				.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("G648", "G654")
				.build();
		protein.flexibility.get("G649").setLibraryRotamers("PHE", "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G650").setLibraryRotamers("ASP").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G651").setLibraryRotamers("GLU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G654").setLibraryRotamers("THR").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A155", "A194")
				.build();
		ligand.flexibility.get("A156").setLibraryRotamers("PHE").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A172").setLibraryRotamers("LYS").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A192").setLibraryRotamers("ILE").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A193").setLibraryRotamers("THR").addWildTypeRotamers().setContinuous();

		// make the conf spaces ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
				.addStrand(ligand)
				.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
				.addStrands(protein, ligand)
				.build();
		double epsilon = 0.01;

		Result result = runKStar(confSpaces, epsilon, null);

	}

	private static class ConfSpaces {
		public ForcefieldParams ffparams;
		public SimpleConfSpace protein;
		public SimpleConfSpace ligand;
		public SimpleConfSpace complex;
	}

	private static class Result {
		public KStar kstar;
		public List<KStar.ScoredSequence> scores;
	}

	public static Result runKStar(ConfSpaces confSpaces, double epsilon, String confDBPattern) {

		AtomicReference<Result> resultRef = new AtomicReference<>(null);

		Parallelism parallelism = Parallelism.makeCpu(4);
		//Parallelism parallelism = Parallelism.make(4, 1, 8);

		// how should we compute energies of molecules?
		new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
			.setParallelism(parallelism)
			.use((ecalc) -> {

				// how should we define energies of conformations?
				KStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
					return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
						.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
							.build()
							.calcReferenceEnergies()
						).build();
				};

				// how should confs be ordered and searched?
				ConfSearchFactory confSearchFactory = (emat, pmat) -> {
					return new ConfAStarTree.Builder(emat, pmat)
						.setTraditional()
						.build();
				};

				KStarScoreWriter.Formatter testFormatter = (KStarScoreWriter.ScoreInfo info) -> {

					Function<PartitionFunction.Result,String> formatPfunc = (result) -> {
						if (result.status == PartitionFunction.Status.Estimated) {
							return String.format("%12e", result.values.qstar.doubleValue());
						}
						return "null";
					};

					return String.format("assertSequence(result, %3d, \"%s\", %-12s, %-12s, %-12s, epsilon); // K* = %s",
						info.sequenceNumber,
						info.sequence.toString(Sequence.Renderer.ResType),
						formatPfunc.apply(info.kstarScore.protein),
						formatPfunc.apply(info.kstarScore.ligand),
						formatPfunc.apply(info.kstarScore.complex),
						info.kstarScore.toString()
					);
				};

				// run K*
				Result result = new Result();
				KStar.Settings settings = new KStar.Settings.Builder()
					.setEpsilon(epsilon)
					.setStabilityThreshold(null)
					.addScoreConsoleWriter()
					.setConfDBPattern(confDBPattern)
					//.setShowPfuncProgress(true)
					.build();
				result.kstar = new KStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, ecalc, confEcalcFactory, confSearchFactory, settings);
				result.scores = result.kstar.run();

				// pass back the ref
				resultRef.set(result);
			});

		return resultRef.get();
	}
}
