package edu.duke.cs.osprey.kstar;

import static edu.duke.cs.osprey.TestBase.fileForWriting;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.KStar.ConfSearchFactory;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.junit.Test;

import java.util.HashSet;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Function;

public class TestKStar {

	public static class ConfSpaces {
		public ForcefieldParams ffparams;
		public SimpleConfSpace protein;
		public SimpleConfSpace ligand;
		public SimpleConfSpace complex;
	}

	public static class Result {
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
					.addScoreConsoleWriter(testFormatter)
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

	public static ConfSpaces make2RL0() {

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
		protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "ASP").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType, "SER", "ASN", "GLN").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType, "ASP", "GLU").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "PHE", "TYR").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType, "SER", "ASN").addWildTypeRotamers().setContinuous();

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

		return confSpaces;
	}

	@Test
	public void test2RL0() {

		double epsilon = 0.95;
		Result result = runKStar(make2RL0(), epsilon, null);
		assert2RL0(result, epsilon);
	}

	private static void assert2RL0(Result result, double epsilon) {
		// check the results (values collected with e = 0.1 and 64 digits precision)
		// NOTE: these values don't match the ones in the TestKSImplLinear test because the conf spaces are slightly different
		// also, the new K* code has been updated to be more precise
		assertSequence(result,   0, "PHE ASP GLU THR PHE LYS ILE THR", 4.300422e+04, 4.347270e+30, 4.201039e+50, epsilon); // K* = 15.351629 in [15.312533,15.396058] (log10)
		assertSequence(result,   1, "PHE ASP GLU THR PHE LYS ILE SER", 4.300422e+04, 1.076556e+30, 4.045744e+50, epsilon); // K* = 15.941451 in [15.878562,15.986656] (log10)
		assertSequence(result,   2, "PHE ASP GLU THR PHE LYS ILE ASN", 4.300422e+04, 4.650623e+29, 1.854792e+49, epsilon); // K* = 14.967273 in [14.920727,15.011237] (log10)
		assertSequence(result,   3, "PHE ASP GLU THR PHE LYS ALA THR", 4.300422e+04, 1.545055e+27, 7.003938e+45, epsilon); // K* = 14.022887 in [13.984844,14.066296] (log10)
		assertSequence(result,   4, "PHE ASP GLU THR PHE LYS VAL THR", 4.300422e+04, 5.694044e+28, 9.854022e+47, epsilon); // K* = 14.604682 in [14.558265,14.649295] (log10)
		assertSequence(result,   5, "PHE ASP GLU THR PHE LYS LEU THR", 4.300422e+04, 3.683508e-11, 3.644143e+08, epsilon); // K* = 14.361823 in [14.324171,14.405292] (log10)
		assertSequence(result,   6, "PHE ASP GLU THR PHE LYS PHE THR", 4.300422e+04, 2.820863e+24, null        , epsilon); // K* = none      in [-Infinity,-Infinity] (log10)
		assertSequence(result,   7, "PHE ASP GLU THR PHE LYS TYR THR", 4.300422e+04, 1.418587e+26, null        , epsilon); // K* = none      in [-Infinity,-Infinity] (log10)
		assertSequence(result,   8, "PHE ASP GLU THR PHE ASP ILE THR", 4.300422e+04, 4.294128e+20, 1.252820e+36, epsilon); // K* = 10.831503 in [10.802450,10.873479] (log10)
		assertSequence(result,   9, "PHE ASP GLU THR PHE GLU ILE THR", 4.300422e+04, 4.831904e+20, 2.273475e+35, epsilon); // K* = 10.039061 in [10.009341,10.081194] (log10)
		assertSequence(result,  10, "PHE ASP GLU THR TYR LYS ILE THR", 4.300422e+04, 4.583429e+30, 2.294671e+50, epsilon); // K* = 15.066019 in [15.026963,15.110683] (log10)
		assertSequence(result,  11, "PHE ASP GLU THR ALA LYS ILE THR", 4.300422e+04, 3.310340e+28, 2.171286e+47, epsilon); // K* = 14.183333 in [14.156555,14.226141] (log10)
		assertSequence(result,  12, "PHE ASP GLU THR VAL LYS ILE THR", 4.300422e+04, 9.004068e+29, 1.866542e+49, epsilon); // K* = 14.683088 in [14.652599,14.725905] (log10)
		assertSequence(result,  13, "PHE ASP GLU THR ILE LYS ILE THR", 4.300422e+04, 3.398648e+30, 1.598348e+50, epsilon); // K* = 15.038854 in [15.002827,15.082763] (log10)
		assertSequence(result,  14, "PHE ASP GLU THR LEU LYS ILE THR", 4.300422e+04, 5.296285e+27, 1.045234e+47, epsilon); // K* = 14.661731 in [14.616778,14.705997] (log10)
		assertSequence(result,  15, "PHE ASP GLU SER PHE LYS ILE THR", 3.153051e+06, 4.347270e+30, 1.477525e+53, epsilon); // K* = 16.032587 in [15.970427,16.078237] (log10)
		assertSequence(result,  16, "PHE ASP GLU ASN PHE LYS ILE THR", 1.484782e+06, 4.347270e+30, 9.591789e+52, epsilon); // K* = 16.172020 in [16.114119,16.217413] (log10)
		assertSequence(result,  17, "PHE ASP GLU GLN PHE LYS ILE THR", 2.531411e+06, 4.347270e+30, 3.346589e+53, epsilon); // K* = 16.483023 in [16.424659,16.528464] (log10)
		assertSequence(result,  18, "PHE ASP ASP THR PHE LYS ILE THR", 1.216972e+01, 4.347270e+30, 1.469949e+45, epsilon); // K* = 13.443805 in [13.413093,13.487750] (log10)
		assertSequence(result,  19, "PHE GLU GLU THR PHE LYS ILE THR", 1.986991e+05, 4.347270e+30, 1.097189e+50, epsilon); // K* = 14.103869 in [14.056796,14.148018] (log10)
		assertSequence(result,  20, "TYR ASP GLU THR PHE LYS ILE THR", 1.666243e+04, 4.347270e+30, 2.814673e+46, epsilon); // K* = 11.589473 in [11.550098,11.633104] (log10)
		assertSequence(result,  21, "ALA ASP GLU THR PHE LYS ILE THR", 6.100779e+02, 4.347270e+30, 1.671418e+45, epsilon); // K* = 11.799483 in [11.777675,11.841368] (log10)
		assertSequence(result,  22, "VAL ASP GLU THR PHE LYS ILE THR", 1.271497e+02, 4.347270e+30, 2.380877e+45, epsilon); // K* = 12.634205 in [12.613207,12.676919] (log10)
		assertSequence(result,  23, "ILE ASP GLU THR PHE LYS ILE THR", 5.942890e+02, 4.347270e+30, 2.012605e+46, epsilon); // K* = 12.891544 in [12.844167,12.936569] (log10)
		assertSequence(result,  24, "LEU ASP GLU THR PHE LYS ILE THR", 4.614233e+00, 4.347270e+30, 4.735376e+43, epsilon); // K* = 12.373038 in [12.339795,12.417250] (log10)
	}

	public static ConfSpaces make1GUA11() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
			.addMoleculeForWildTypeRotamers(mol)
			.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("1", "180")
			.build();
		protein.flexibility.get("21").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("24").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("25").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("27").setLibraryRotamers(Strand.WildType, "HID").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("29").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("40").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("181", "215")
			.build();
		ligand.flexibility.get("209").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("213").setLibraryRotamers(Strand.WildType, "HID", "HIE", "LYS", "ARG").addWildTypeRotamers().setContinuous();

		// make the complex conf space ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
			.addStrand(protein)
			.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
			.addStrands(protein, ligand)
			.build();

		return confSpaces;
	}

	@Test
	public void test1GUA11() {

		double epsilon = 0.999999;
		Result result = runKStar(make1GUA11(), epsilon, null);

		// check the results (values collected with e = 0.1 and 64 digits precision)
		assertSequence(result,   0, "ILE ILE GLN HIE VAL TYR LYS VAL", 1.186071e+42, 2.840001e+07, 1.119884e+66, epsilon); // K* = 16.521744 in [16.463680,16.563832] (log10)
		assertSequence(result,   1, "ILE ILE GLN HIE VAL TYR LYS HID", 1.186071e+42, 5.575412e+07, 3.345731e+66, epsilon); // K* = 16.704103 in [16.647717,16.747742] (log10)
		assertSequence(result,   2, "ILE ILE GLN HIE VAL TYR LYS HIE", 1.186071e+42, 5.938851e+06, 5.542993e+65, epsilon); // K* = 16.895931 in [16.826906,16.938784] (log10)
		assertSequence(result,   3, "ILE ILE GLN HIE VAL TYR LYS LYS", 1.186071e+42, 6.402058e+04, 3.315165e+63, epsilon); // K* = 16.640075 in [16.563032,16.685734] (log10)
		assertSequence(result,   4, "ILE ILE GLN HIE VAL TYR LYS ARG", 1.186071e+42, 1.157637e+05, 5.375731e+64, epsilon); // K* = 17.592754 in [17.514598,17.638466] (log10)
		assertSequence(result,   5, "ILE ILE GLN HID VAL TYR LYS VAL", 9.749716e+41, 2.840001e+07, 2.677894e+66, epsilon); // K* = 16.985483 in [16.927677,17.026890] (log10)
	}

	@Test
	public void test2RL0WithConfDB() {

		final double epsilon = 0.95;
		final String confdbPattern = "kstar.*.conf.db";
		final ConfSpaces confSpaces = make2RL0();

		fileForWriting("kstar.protein.conf.db", (proteinDBFile) -> {
			fileForWriting("kstar.ligand.conf.db", (ligandDBFile) -> {
				fileForWriting("kstar.complex.conf.db", (complexDBFile) -> {

					// run with empty dbs
					Stopwatch sw = new Stopwatch().start();
					Result result = runKStar(confSpaces, epsilon, confdbPattern);
					assert2RL0(result, epsilon);
					System.out.println(sw.getTime(2));

					// the dbs should have stuff in them

					new ConfDB(confSpaces.protein, proteinDBFile).use((confdb) -> {
						HashSet<Sequence> sequences = new HashSet<>(result.kstar.protein.sequences);
						assertThat(confdb.getNumSequences(), is((long)sequences.size()));
						for (Sequence sequence : confdb.getSequences()) {
							assertThat(confdb.getSequence(sequence).size(), greaterThan(0L));
						}
					});

					new ConfDB(confSpaces.ligand, ligandDBFile).use((confdb) -> {
						HashSet<Sequence> sequences = new HashSet<>(result.kstar.ligand.sequences);
						assertThat(confdb.getNumSequences(), is((long)sequences.size()));
						for (Sequence sequence : confdb.getSequences()) {
							assertThat(confdb.getSequence(sequence).size(), greaterThan(0L));
						}
					});

					new ConfDB(confSpaces.complex, complexDBFile).use((confdb) -> {
						List<Sequence> sequences = result.kstar.complex.sequences;
						assertThat(confdb.getNumSequences(), is((long)sequences.size()));
						for (Sequence sequence : confdb.getSequences()) {
							assertThat(confdb.getSequence(sequence).size(), greaterThan(0L));
						}
					});

					assertThat(proteinDBFile.exists(), is(true));
					assertThat(ligandDBFile.exists(), is(true));
					assertThat(complexDBFile.exists(), is(true));

					// run again with full dbs
					sw = new Stopwatch().start();
					Result result2 = runKStar(confSpaces, epsilon, confdbPattern);
					assert2RL0(result2, epsilon);
					System.out.println(sw.getTime(2));
				});
			});
		});
	}

	public static void assertSequence(Result result, int sequenceIndex, String sequence, Double proteinQStar, Double ligandQStar, Double complexQStar, double epsilon) {

		KStar.ScoredSequence scoredSequence = result.scores.get(sequenceIndex);

		// check the sequence
		assertThat(scoredSequence.sequence.toString(Sequence.Renderer.ResType), is(sequence));

		// check q* values and epsilon
		assertResult(scoredSequence.score.protein, proteinQStar, epsilon);
		assertResult(scoredSequence.score.ligand, ligandQStar, epsilon);
		assertResult(scoredSequence.score.complex, complexQStar, epsilon);
	}

	public static void assertResult(PartitionFunction.Result result, Double qstar, double epsilon) {
		if (qstar != null) {
			assertThat(result.status, is(PartitionFunction.Status.Estimated));
			assertThat(result.values.qstar.doubleValue(), greaterThanOrEqualTo(qstar*(1.0 - epsilon)));
			assertThat(result.values.getEffectiveEpsilon(), lessThanOrEqualTo(epsilon));
		} else {
			assertThat(result.status, is(not(PartitionFunction.Status.Estimated)));
		}
	}
}
