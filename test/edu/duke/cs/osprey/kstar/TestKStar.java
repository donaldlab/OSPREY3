package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
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
import org.junit.Test;

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
		public List<KStarScore> scores;
	}

	public static Result runKStar(ConfSpaces confSpaces, double epsilon) {

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

					return String.format("assertSequence(result, %3d, \"%s\", %-12s, %-12s, %-12s, epsilon); // K*(log10) = %s",
						info.sequenceNumber,
						info.sequence,
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
					.addScoreConsoleWriter(testFormatter)
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
			.setResidues(648, 654)
			.build();
		protein.flexibility.get(649).setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get(650).setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get(651).setLibraryRotamers(Strand.WildType, "ASP").addWildTypeRotamers().setContinuous();
		protein.flexibility.get(654).setLibraryRotamers(Strand.WildType, "SER", "ASN", "GLN").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues(155, 194)
			.build();
		ligand.flexibility.get(156).setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get(172).setLibraryRotamers(Strand.WildType, "ASP", "GLU").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get(192).setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "PHE", "TYR").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get(193).setLibraryRotamers(Strand.WildType, "SER", "ASN").addWildTypeRotamers().setContinuous();

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
		Result result = runKStar(make2RL0(), epsilon);

		// check the results (values collected with e = 0.1)
		// NOTE: these values don't match the ones in the TestKSImplLinear test because the conf spaces are slightly different
		assertSequence(result,   0, "PHE ASP GLU THR PHE LYS ILE THR", 4.315645e+04, 4.347270e+30, 4.201039e+50, epsilon); // K*(log10) = 15.350094 in [15.315503,15.394524]
		assertSequence(result,   1, "PHE ASP GLU THR PHE LYS ILE SER", 4.315645e+04, 1.076556e+30, 4.045744e+50, epsilon); // K*(log10) = 15.939916 in [15.881531,15.985122]
		assertSequence(result,   2, "PHE ASP GLU THR PHE LYS ILE ASN", 4.315645e+04, 4.650623e+29, 1.854801e+49, epsilon); // K*(log10) = 14.965741 in [14.923698,15.009355]
		assertSequence(result,   3, "PHE ASP GLU THR PHE LYS ALA THR", 4.315645e+04, 1.545055e+27, 7.003938e+45, epsilon); // K*(log10) = 14.021353 in [13.987814,14.064762]
		assertSequence(result,   4, "PHE ASP GLU THR PHE LYS VAL THR", 4.315645e+04, 5.694044e+28, 9.854022e+47, epsilon); // K*(log10) = 14.603147 in [14.561234,14.647761]
		assertSequence(result,   5, "PHE ASP GLU THR PHE LYS LEU THR", 4.315645e+04, null        , 3.644115e+08, epsilon); // K*(log10) = none in [none,none]
		assertSequence(result,   6, "PHE ASP GLU THR PHE LYS PHE THR", 4.315645e+04, 2.820863e+24, null        , epsilon); // K*(log10) = none in [-Infinity,-Infinity]
		assertSequence(result,   7, "PHE ASP GLU THR PHE LYS TYR THR", 4.315645e+04, 1.418587e+26, null        , epsilon); // K*(log10) = none in [-Infinity,-Infinity]
		assertSequence(result,   8, "PHE ASP GLU THR PHE ASP ILE THR", 4.315645e+04, 4.289012e+20, 1.252934e+36, epsilon); // K*(log10) = 10.830525 in [10.805317,10.871671]
		assertSequence(result,   9, "PHE ASP GLU THR PHE GLU ILE THR", 4.315645e+04, 4.831904e+20, 2.273475e+35, epsilon); // K*(log10) = 10.037526 in [10.012310,10.079659]
		assertSequence(result,  10, "PHE ASP GLU THR TYR LYS ILE THR", 4.315645e+04, 4.583429e+30, 2.294685e+50, epsilon); // K*(log10) = 15.064487 in [15.029935,15.108900]
		assertSequence(result,  11, "PHE ASP GLU THR ALA LYS ILE THR", 4.315645e+04, 3.305184e+28, 2.171286e+47, epsilon); // K*(log10) = 14.182476 in [14.159251,14.225284]
		assertSequence(result,  12, "PHE ASP GLU THR VAL LYS ILE THR", 4.315645e+04, 9.004068e+29, 1.866542e+49, epsilon); // K*(log10) = 14.681553 in [14.655568,14.724370]
		assertSequence(result,  13, "PHE ASP GLU THR ILE LYS ILE THR", 4.315645e+04, 3.398648e+30, 1.598348e+50, epsilon); // K*(log10) = 15.037320 in [15.005796,15.081229]
		assertSequence(result,  14, "PHE ASP GLU THR LEU LYS ILE THR", 4.315645e+04, 5.296285e+27, 1.045234e+47, epsilon); // K*(log10) = 14.660196 in [14.619748,14.704462]
		assertSequence(result,  15, "PHE ASP GLU SER PHE LYS ILE THR", 3.153051e+06, 4.347270e+30, 1.477525e+53, epsilon); // K*(log10) = 16.032587 in [15.970427,16.078237]
		assertSequence(result,  16, "PHE ASP GLU ASN PHE LYS ILE THR", 1.487639e+06, 4.347270e+30, 9.591789e+52, epsilon); // K*(log10) = 16.171185 in [16.114634,16.216578]
		assertSequence(result,  17, "PHE ASP GLU GLN PHE LYS ILE THR", 2.531411e+06, 4.347270e+30, 3.346589e+53, epsilon); // K*(log10) = 16.483023 in [16.424659,16.528464]
		assertSequence(result,  18, "PHE ASP ASP THR PHE LYS ILE THR", 1.211611e+01, 4.347270e+30, 1.469961e+45, epsilon); // K*(log10) = 13.445726 in [13.412206,13.489261]
		assertSequence(result,  19, "PHE GLU GLU THR PHE LYS ILE THR", 1.986991e+05, 4.347270e+30, 1.097189e+50, epsilon); // K*(log10) = 14.103869 in [14.056796,14.148018]
		assertSequence(result,  20, "TYR ASP GLU THR PHE LYS ILE THR", 1.666243e+04, 4.347270e+30, 2.814673e+46, epsilon); // K*(log10) = 11.589473 in [11.550098,11.633104]
		assertSequence(result,  21, "ALA ASP GLU THR PHE LYS ILE THR", 6.113463e+02, 4.347270e+30, 1.671418e+45, epsilon); // K*(log10) = 11.798581 in [11.778039,11.840466]
		assertSequence(result,  22, "VAL ASP GLU THR PHE LYS ILE THR", 1.269943e+02, 4.347270e+30, 2.380877e+45, epsilon); // K*(log10) = 12.634736 in [12.612457,12.677450]
		assertSequence(result,  23, "ILE ASP GLU THR PHE LYS ILE THR", 5.942890e+02, 4.347270e+30, 2.012585e+46, epsilon); // K*(log10) = 12.891540 in [12.844163,12.936699]
		assertSequence(result,  24, "LEU ASP GLU THR PHE LYS ILE THR", 4.614195e+00, 4.347270e+30, 4.735376e+43, epsilon); // K*(log10) = 12.373042 in [12.336269,12.417254]
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
			.setResidues(1, 180)
			.build();
		protein.flexibility.get(21).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get(24).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get(25).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get(27).setLibraryRotamers(Strand.WildType, "HID").addWildTypeRotamers().setContinuous();
		protein.flexibility.get(29).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get(40).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues(181, 215)
			.build();
		ligand.flexibility.get(209).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get(213).setLibraryRotamers(Strand.WildType, "HID", "HIE", "LYS", "ARG").addWildTypeRotamers().setContinuous();

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
	public void test1GUA_11() {

		double epsilon = 0.99;
		Result result = runKStar(make1GUA11(), epsilon);

		// check the results (values collected with e = 0.1)
		assertSequence(result,   0, "ILE ILE GLN HIE VAL TYR LYS VAL", 1.1838e+42, 2.7098e+7, 1.1195e+66, epsilon);
		assertSequence(result,   1, "ILE ILE GLN HIE VAL TYR LYS HID", 1.1838e+42, 5.5334e+7, 3.3455e+66, epsilon);
		assertSequence(result,   2, "ILE ILE GLN HIE VAL TYR LYS HIE", 1.1838e+42, 5.8485e+6, 5.5426e+65, epsilon);
		assertSequence(result,   3, "ILE ILE GLN HIE VAL TYR LYS LYS", 1.1838e+42, 6.3856e+4, 3.3162e+63, epsilon);
		assertSequence(result,   4, "ILE ILE GLN HIE VAL TYR LYS ARG", 1.1838e+42, 1.1527e+5, 5.3772e+64, epsilon);
		assertSequence(result,   5, "ILE ILE GLN HID VAL TYR LYS VAL", 9.7448e+41, 2.7098e+7, 2.6775e+66, epsilon);
	}

	public static void assertSequence(Result result, int sequenceIndex, String sequence, Double proteinQStar, Double ligandQStar, Double complexQStar, double epsilon) {

		KStarScore score = result.scores.get(sequenceIndex);

		// check the sequence
		assertThat(result.kstar.complex.sequences.get(sequenceIndex), is(new KStar.Sequence(sequence)));

		// check q* values and epsilon
		assertResult(score.protein, proteinQStar, epsilon);
		assertResult(score.ligand, ligandQStar, epsilon);
		assertResult(score.complex, complexQStar, epsilon);
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
