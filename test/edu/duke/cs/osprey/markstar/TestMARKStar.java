package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.TestBBKStar;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.MARKStar.ConfSearchFactory;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.junit.Test;

import java.io.File;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.List;

import static edu.duke.cs.osprey.kstar.TestBBKStar.runBBKStar;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.assertThat;

//import edu.duke.cs.osprey.kstar.KStar.ConfSearchFactory;

public class TestMARKStar {

	public static final int NUM_CPUs = 4;
	public static boolean REUDCE_MINIMIZATIONS = false;
	public static final EnergyPartition ENERGY_PARTITION = EnergyPartition.Traditional;

	public static class ConfSpaces {
		public ForcefieldParams ffparams;
		public SimpleConfSpace protein;
		public SimpleConfSpace ligand;
		public SimpleConfSpace complex;
	}

	public static class Result {
		public MARKStar markstar;
		public List<MARKStar.ScoredSequence> scores;
	}

	@Test
	public void test3K3Q () {

		ConfSpaces confSpaces = makeConfSpaces(
			"examples/python.KStar/3k3q_prepped.pdb",
				new String[]{"C253", "C417"},
				new String[]{"C412", "C404", "C409", "C391", "C390", "C384", "C392", "C383", "C397", "C378", "C376", "C413"},
				new String[]{"B2","B250"},
				new String[]{"B193", "B216", "B214"}
		);
		final double epsilon = 0.9999;
		String kstartime = "(not run)";
		boolean runkstar = true;
		Stopwatch runtime = new Stopwatch().start();
		if(runkstar) {
			List<KStar.ScoredSequence> seqs = runKStar(confSpaces, epsilon);
			runtime.stop();
			kstartime = runtime.getTime(2);
			runtime.reset();
			runtime.start();
		}
		Result result = runMARKStar(confSpaces, epsilon);
		runtime.stop();
		String markstartime = runtime.getTime(2);
		for(MARKStar.ScoredSequence seq: result.scores)
			printMARKStarComputationStats(seq);
		System.out.println("MARK* time: "+markstartime+", K* time: "+kstartime);
	}

	private ConfSpaces makeConfSpaces(String pdb, String[] proteinDef, String[] proteinFlex, String[] ligandDef,
								String[] ligandFlex) {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readFile(pdb);

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
			.addMoleculeForWildTypeRotamers(mol)
			.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues(proteinDef[0], proteinDef[1])
			.build();
		for(String resName: proteinFlex)
            protein.flexibility.get(resName).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues(ligandDef[0], ligandDef[1])
			.build();
		for(String resName: ligandFlex)
			ligand.flexibility.get(resName).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

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
	public void test4KT6 () {
		ConfSpaces confSpaces = make4KT6();
		final double epsilon = 0.99;
		String kstartime = "(not run)";
		boolean runkstar = true;
		Stopwatch runtime = new Stopwatch().start();
		if(runkstar) {
			List<KStar.ScoredSequence> seqs = runKStar(confSpaces, epsilon);
			runtime.stop();
			kstartime = runtime.getTime(2);
			runtime.reset();
			runtime.start();
		}
		Result result = runMARKStar(confSpaces, epsilon);
		runtime.stop();
		String markstartime = runtime.getTime(2);
		for(MARKStar.ScoredSequence seq: result.scores)
			printMARKStarComputationStats(seq);
		System.out.println("MARK* time: "+markstartime+", K* time: "+kstartime);
	}

	private ConfSpaces make4KT6() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readFile("examples/python.KStar/4kt6_prepped.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
			.addMoleculeForWildTypeRotamers(mol)
			.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("C193", "C446")
			.build();
		protein.flexibility.get("C290").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("D1", "D161")
			.build();
		ligand.flexibility.get("D148").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("D149").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("D151").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("D152").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("D153").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("D155").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("D156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

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
	public void test2XXM() {
		ConfSpaces confSpaces = make2XXM();
		final double epsilon = 0.99;
		String kstartime = "(not run)";
		boolean runkstar = true;
		Stopwatch runtime = new Stopwatch().start();
		List<KStar.ScoredSequence> kStarSeqs = null;
		if(runkstar) {
			kStarSeqs = runKStar(confSpaces, epsilon);
			runtime.stop();
			kstartime = runtime.getTime(2);
			runtime.reset();
			runtime.start();
		}
		Result result = runMARKStar(confSpaces, epsilon);
		runtime.stop();
		String markstartime = runtime.getTime(2);
		System.out.println("MARK* time: "+markstartime+", K* time: "+kstartime);
		for(MARKStar.ScoredSequence seq: result.scores)
			printMARKStarComputationStats(seq);
		if(runkstar)
			for(KStar.ScoredSequence seq: kStarSeqs)
				printKStarComputationStats(seq);
	}

	private ConfSpaces make2XXM() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readFile("examples/python.KStar/2xxm_prepped.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.addMoleculeForWildTypeRotamers(mol)
				.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A146", "A218")
				.build();
		protein.flexibility.get("A177").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A178").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A179").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A180").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("A181").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("B4", "B113")
				.build();
		ligand.flexibility.get("B58").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("B60").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("B61").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("B64").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

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
	public void test1A0RBBKStar() {
		TestKStar.ConfSpaces confSpaces = make1A0RBBKStar();
		int numSequences = 2;
		double epsilon = 0.999;
		String kstartime = "(Not run)";
		Stopwatch watch = new Stopwatch();
		boolean runkstar = true;
		if(runkstar) {
			watch.start();
			runBBKStar(confSpaces, numSequences, epsilon, null, 1, false);
			watch.stop();
			kstartime = watch.getTime(2);
			watch.reset();
		}
		watch.start();
		runBBKStar(confSpaces, numSequences, epsilon, null, 1, true);
		watch.stop();
		String bbkstartime = watch.getTime(2);
		System.out.println("MARK*: "+bbkstartime+", Traditional K*: "+kstartime);
	}

	private static TestKStar.ConfSpaces make1A0RBBKStar() {

		TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readFile("examples/python.KStar/1a0r_prepped.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("P13", "P230")
				.build();
		protein.flexibility.get("P219").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("P220").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("P222").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("P223").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "ILE", "PHE", "TYR",
				"TRP", "CYS", "MET", "SER", "THR", "LYS", "ARG", "HIP", "HIE", "HID", "ASP", "GLU", "ASN", "GLN",
				"GLY").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("P224").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();


		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("B2", "B340")
				.build();
		ligand.flexibility.get("B46").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("B47").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

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
	public void test1GUA11MARKVsTraditional() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make1GUA11();
		final double epsilon = 0.999;
		final int numSequences = 6;


		String traditionalTime = "(Not run)";
		boolean runkstar = false;
		Stopwatch timer = new Stopwatch().start();
		if(runkstar) {
			TestBBKStar.Results results = runBBKStar(confSpaces, numSequences, epsilon, null, 1, false);
			timer.stop();
			traditionalTime = timer.getTime(2);
			timer.reset();
			timer.start();
		}
		runBBKStar(confSpaces, numSequences, epsilon, null, 1, true);
		String MARKStarTime = timer.getTime(2);
		timer.stop();

		//assert2RL0(results, numSequences);
		System.out.println("Traditional time: "+traditionalTime);
		System.out.println("MARK* time: "+MARKStarTime);
	}


	@Test
	public void timeMARKStarVsTraditional() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();
		final double epsilon = 0.68;
		final int numSequences = 10;
		Stopwatch timer = new Stopwatch().start();
		TestBBKStar.Results results = runBBKStar(confSpaces, numSequences, epsilon, null, 1, false);
		timer.stop();
		String traditionalTime = timer.getTime(2);
		timer.reset();
		timer.start();
		results = runBBKStar(confSpaces, numSequences, epsilon, null, 1, true);
		String MARKStarTime = timer.getTime(2);
		timer.stop();

		//assert2RL0(results, numSequences);
		System.out.println("Traditional time: "+traditionalTime);
		System.out.println("MARK* time: "+MARKStarTime);
	}

	@Test
	public void testMARKStarVsKStar() {
		int numFlex = 10;
		double epsilon = 0.68;
		compareMARKStarAndKStar(numFlex, epsilon);
	}

	@Test
	public void testMARKStarTinyEpsilon() {
		printMARKStarComputationStats(runMARKStar(3, 0.9).get(0));

	}

	@Test
	public void compareReducedMinimizationsVsNormal() {
		int numFlex = 14;
		double epsilon = 0.68;
		REUDCE_MINIMIZATIONS = false;
		System.out.println("Trying without reduced minimizations...");
		Stopwatch runTime = new Stopwatch().start();
		runMARKStar(numFlex, epsilon);
		String withoutRMTime= runTime.getTime(2);
		runTime.stop();
		runTime.reset();
		runTime.start();
		REUDCE_MINIMIZATIONS = true;
		System.out.println("Retrying with reduced minimizations...");
		runMARKStar(numFlex, epsilon);
		runTime.stop();
		System.out.println("Without Reduced Minimization time: "+withoutRMTime);
		System.out.println("Reduced minimization time: "+runTime.getTime(2));
	}

	@Test
	public void test1GUASmallUpTo()
	{
		int maxNumFlex = 8;
		double epsilon = 0.68;
		for(int i = 1; i < maxNumFlex; i++)
			compareMARKStarAndKStar(i,0.68);
	}

	private void compareMARKStarAndKStar(int numFlex, double epsilon) {
		Stopwatch runTime = new Stopwatch().start();
		String kstartime = "(not run)";
		List<KStar.ScoredSequence> kStarSeqs = null;
		boolean runkstar = true;
		if(runkstar) {
			kStarSeqs = runKStarComparison(numFlex, epsilon);
			runTime.stop();
			kstartime = runTime.getTime(2);
			runTime.reset();
			runTime.start();
		}
		List<MARKStar.ScoredSequence> markStarSeqs = runMARKStar(numFlex, epsilon);
		runTime.stop();
		for(MARKStar.ScoredSequence seq: markStarSeqs)
			printMARKStarComputationStats(seq);
		if(runkstar)
			for(KStar.ScoredSequence seq: kStarSeqs)
				printKStarComputationStats(seq);
		System.out.println("K* time: "+kstartime);
		System.out.println("MARK* time: "+runTime.getTime(2));
	}


	private void printMARKStarComputationStats(MARKStar.ScoredSequence result) {
		int totalConfsEnergied = result.score.complex.numConfs + result.score.protein.numConfs + result.score.ligand.numConfs;
		int totalConfsLooked = result.score.complex.getNumConfsLooked()+ result.score.protein.getNumConfsLooked()+ result.score.ligand.getNumConfsLooked();
		BigInteger totalConfSpaceSize = new BigInteger(result.score.complex.totalNumConfs)
				.add(new BigInteger(result.score.protein.totalNumConfs)).add(new BigInteger(result.score.ligand.totalNumConfs));
		System.out.println("MARK* Stats: "+String.format("score:%12e in [%12e,%12e] (log10), confs looked at:%4d, confs minimized:%4d\n total confSize:%4s",MathTools.log10p1(result.score.score), MathTools.log10p1(result.score.lowerBound),
				MathTools.log10p1(result.score.upperBound),totalConfsLooked,totalConfsEnergied, totalConfSpaceSize.toString()));
		System.out.println("Above stats for sequence: "+result.sequence);

	}

	private void printKStarComputationStats(KStar.ScoredSequence result) {
		int totalConfsEnergied = result.score.complex.numConfs + result.score.protein.numConfs + result.score.ligand.numConfs;
		int totalConfsLooked = result.score.complex.getNumConfsLooked()+ result.score.protein.getNumConfsLooked()+ result.score.ligand.getNumConfsLooked();
		BigInteger totalConfSpaceSize = new BigInteger(result.score.complex.totalNumConfs)
				.add(new BigInteger(result.score.protein.totalNumConfs)).add(new BigInteger(result.score.ligand.totalNumConfs));
		System.out.println("K* Stats: "+String.format("score:%12e in [%12e,%12e] (log10), confs looked at:%4d, confs minimized:%4d\ntotal confSize:%4s",MathTools.log10p1(result.score.score), MathTools.log10p1(result.score.lowerBound),
				MathTools.log10p1(result.score.upperBound),totalConfsLooked,totalConfsEnergied, totalConfSpaceSize.toString()));
		System.out.println("Above stats for sequence: "+result.sequence);
	}

	@Test
	public void testMARKStar2RL0(){

	}

	@Test
	public void KStarComparison() {
		List<KStar.ScoredSequence> results = runKStarComparison(5,0.68);
		for (int index = 0; index < results.size(); index++) {
			int totalConfsEnergied = results.get(index).score.complex.numConfs + results.get(index).score.protein.numConfs + results.get(index).score.ligand.numConfs;
			int totalConfsLooked = results.get(index).score.complex.getNumConfsLooked()+ results.get(index).score.protein.getNumConfsLooked()+ results.get(index).score.ligand.getNumConfsLooked();
			System.out.println(String.format("score:%12e in [%12e,%12e], confs looked at:%4d, confs minimized:%4d",results.get(index).score.score, results.get(index).score.lowerBound,
					results.get(index).score.upperBound,totalConfsLooked,totalConfsEnergied));
		}
	}



	@Test
	public void testMARKStar() {
		List<MARKStar.ScoredSequence> results = runMARKStar(1, 0.01);
		for (int index = 0; index < results.size(); index++) {
			int totalConfsEnergied = results.get(index).score.complex.numConfs + results.get(index).score.protein.numConfs + results.get(index).score.ligand.numConfs;
			int totalConfsLooked = results.get(index).score.complex.getNumConfsLooked()+ results.get(index).score.protein.getNumConfsLooked()+ results.get(index).score.ligand.getNumConfsLooked();
			System.out.println(String.format("score:%12e in [%12e,%12e], confs looked at:%4d, confs minimized:%4d",results.get(index).score.score, results.get(index).score.lowerBound,
					results.get(index).score.upperBound,totalConfsLooked,totalConfsEnergied));
		}
	}

	private static List<MARKStar.ScoredSequence> runMARKStar(int numFlex, double epsilon) {
		//ConfSpaces confSpaces = make1GUASmallCATS(numFlex);
		//ConfSpaces confSpaces = make1GUASmallDEEP(numFlex);
		ConfSpaces confSpaces = make1GUASmall(numFlex);
		Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

		// Define the minimizing energy calculator
		EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
				.setParallelism(parallelism)
				.setIsMinimizing(true)
				.build();
		// Define the rigid energy calculator
		EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
				.setIsMinimizing(false)
				.build();
		// how should we define energies of conformations?
		MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
			return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, minimizingEcalc)
							.build()
							.calcReferenceEnergies()
					)
					.setEnergyPartition(ENERGY_PARTITION)
					.build();
		};

		MARKStar.Settings settings = new MARKStar.Settings.Builder()
				.setEpsilon(epsilon)
				.setEnergyMatrixCachePattern("*testmat.emat")
				.setShowPfuncProgress(true)
				.setParallelism(parallelism)
				.setReduceMinimizations(REUDCE_MINIMIZATIONS)
				.build();
		MARKStar run = new MARKStar(confSpaces.protein, confSpaces.ligand,
				confSpaces.complex, rigidEcalc, minimizingEcalc, confEcalcFactory, settings);
		return run.run();
	}

	public static class KstarResult {
		public KStar kstar;
		public List<KStar.ScoredSequence> scores;
	}

	public static List<KStar.ScoredSequence> runKStarComparison(int numFlex, double epsilon) {
		//ConfSpaces confSpaces = make1GUASmallCATS(numFlex);
		//ConfSpaces confSpaces = make1GUASmallDEEP(numFlex);
		ConfSpaces confSpaces = make1GUASmall(numFlex);
		Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

		// Define the minimizing energy calculator
		EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
				.setParallelism(parallelism)
				.build();
		// configure K*
		KStar.Settings settings = new KStar.Settings.Builder()
				.setEpsilon(epsilon)
				.setStabilityThreshold(null)
				.setShowPfuncProgress(true)
				.build();
		KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, settings);
		for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

			// how should we define energies of conformations?
			info.confEcalc = new ConfEnergyCalculator.Builder(info.confSpace, minimizingEcalc)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(info.confSpace, minimizingEcalc)
							.build()
							.calcReferenceEnergies()
					)
					.build();

			// calc energy matrix
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
					.build()
					.calcEnergyMatrix();

			// how should confs be ordered and searched?
			info.confSearchFactory = (rcs) -> {
				ConfAStarTree.Builder builder = new ConfAStarTree.Builder(emat, rcs)
						.setTraditional();
				return builder.build();
			};
		}

		// run K*
		KstarResult result = new KstarResult();
		result.kstar = kstar;
		result.scores = kstar.run();
		return result.scores;
	}

	public static List<KStar.ScoredSequence> runKStar(ConfSpaces confSpaces, double epsilon) {

		Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

		// Define the minimizing energy calculator
		EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
				.setParallelism(parallelism)
				.build();
		// configure K*
		KStar.Settings settings = new KStar.Settings.Builder()
				.setEpsilon(epsilon)
				.setStabilityThreshold(null)
				.setShowPfuncProgress(true)
				.build();
		KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, settings);
		for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

			// how should we define energies of conformations?
			info.confEcalc = new ConfEnergyCalculator.Builder(info.confSpace, minimizingEcalc)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(info.confSpace, minimizingEcalc)
							.build()
							.calcReferenceEnergies()
					)
					.build();

			// calc energy matrix
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
					.setCacheFile(new File(info.type+".kstar.emat"))
					.build()
					.calcEnergyMatrix();

			// how should confs be ordered and searched?
			info.confSearchFactory = (rcs) -> {
				ConfAStarTree.Builder builder = new ConfAStarTree.Builder(emat, rcs)
						.setTraditional();
				return builder.build();
			};
		}

		// run K*
		KstarResult result = new KstarResult();
		result.kstar = kstar;
		result.scores = kstar.run();
		return result.scores;
	}

	public static Result runMARKStar(ConfSpaces confSpaces, double epsilon){

		Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

		// Define the minimizing energy calculator
		EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
				.setParallelism(parallelism)
				.build();
		// Define the rigid energy calculator
		EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
				.setParallelism(parallelism)
				.setIsMinimizing(false)
				.build();
		// how should we define energies of conformations?
		MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
			return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, minimizingEcalc)
							.setCacheFile(new File("test.eref.emat"))
							.build()
							.calcReferenceEnergies()
					)
					.build();
		};

		// how should confs be ordered and searched?
		ConfSearchFactory confSearchFactory = (emat, pmat) -> {
			return new ConfAStarTree.Builder(emat, pmat)
					.setTraditional()
					.build();
		};

		Result result = new Result();

		MARKStar.Settings settings = new MARKStar.Settings.Builder()
				.setEpsilon(epsilon)
				.setEnergyMatrixCachePattern("*testmat.emat")
				.setShowPfuncProgress(true)
				.setParallelism(parallelism)
				.setReduceMinimizations(REUDCE_MINIMIZATIONS)
				.build();

		result.markstar = new MARKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, rigidEcalc, minimizingEcalc, confEcalcFactory, settings);
		result.scores = result.markstar.run();
		return result;
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
		Result result = runMARKStar(make2RL0(), epsilon);

		// check the results (values collected with e = 0.1 and 64 digits precision)
		// NOTE: these values don't match the ones in the TestKSImplLinear test because the conf spaces are slightly different
		// also, the new K* code has been updated to be more precise
		assertSequence(result,   0, "PHE ASP GLU THR PHE LYS ILE THR", 4.300422e+04, 4.347270e+30, 4.201039e+50, epsilon); // K* = 15.351629 in [15.312533,15.396058] (log10)
		assertSequence(result,   1, "PHE ASP GLU THR PHE LYS ILE SER", 4.300422e+04, 1.076556e+30, 4.045744e+50, epsilon); // K* = 15.941451 in [15.878562,15.986656] (log10)
		assertSequence(result,   2, "PHE ASP GLU THR PHE LYS ILE ASN", 4.300422e+04, 4.650623e+29, 1.854792e+49, epsilon); // K* = 14.967273 in [14.920727,15.011237] (log10)
		assertSequence(result,   3, "PHE ASP GLU THR PHE LYS ALA THR", 4.300422e+04, 1.545055e+27, 7.003938e+45, epsilon); // K* = 14.022887 in [13.984844,14.066296] (log10)
		assertSequence(result,   4, "PHE ASP GLU THR PHE LYS VAL THR", 4.300422e+04, 5.694044e+28, 9.854022e+47, epsilon); // K* = 14.604682 in [14.558265,14.649295] (log10)
		//assertSequence(result,   5, "PHE ASP GLU THR PHE LYS LEU THR", 4.300422e+04, 3.683508e-11, 3.644143e+08, epsilon); // K* = 14.361823 in [14.324171,14.405292] (log10)
		//assertSequence(result,   6, "PHE ASP GLU THR PHE LYS PHE THR", 4.300422e+04, 2.820863e+24, null        , epsilon); // K* = none      in [-Infinity,-Infinity] (log10)
		//assertSequence(result,   7, "PHE ASP GLU THR PHE LYS TYR THR", 4.300422e+04, 1.418587e+26, null        , epsilon); // K* = none      in [-Infinity,-Infinity] (log10)
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
	public static ConfSpaces make1GUASmall(int numFlex) {

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
		int start = 21;
		for(int i = start; i < start+numFlex; i++) {
			protein.flexibility.get(i+"").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		}



		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("181", "215")
				.build();
		ligand.flexibility.get("209").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

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
	public static ConfSpaces make1GUASmallDEEP(int numFlex) {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.addMoleculeForWildTypeRotamers(mol)
				.build();

		ArrayList <String> bbflexlist = new ArrayList();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("1", "180")
				.build();
		int start = 21;
		for(int i = start; i < start+numFlex; i++) {
			protein.flexibility.get(i+"").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			bbflexlist.add(i+"");
		}



		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("181", "215")
				.build();
		ligand.flexibility.get("209").setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		bbflexlist.add("209");
		DEEPerSettings deepersettings = new DEEPerSettings(true, "test_deeper.pert", true, "None", false, 2.5,2.5, false, bbflexlist, "/1gua_adj.min.pdb", false, templateLib);
		DEEPerStrandFlex ligand_bbflex = new DEEPerStrandFlex(ligand, deepersettings);

		// make the complex conf space ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
				.addStrand(ligand, ligand_bbflex)
				.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.addStrand(ligand, ligand_bbflex)
				.build();

		return confSpaces;
	}public static ConfSpaces make1GUASmallCATS(int numFlex) {

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
		int start = 21;
		for(int i = start; i < start+numFlex; i++) {
			protein.flexibility.get(i+"").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		}
		CATSStrandFlex bbflex = new CATSStrandFlex(protein, "22", "25");



		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("181", "215")
				.build();
		ligand.flexibility.get("209").setLibraryRotamers(Strand.WildType).addWildTypeRotamers();

		// make the complex conf space ("complex" SimpleConfSpace, har har!)
		confSpaces.protein = new SimpleConfSpace.Builder()
				.addStrand(protein, bbflex)
				.build();
		confSpaces.ligand = new SimpleConfSpace.Builder()
				.addStrand(ligand)
				.build();
		confSpaces.complex = new SimpleConfSpace.Builder()
				.addStrand(protein, bbflex)
				.addStrand(ligand)
				.build();

		return confSpaces;
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
	public void test2RL0Seq5() {
		double epsilon = 0.95;
		Result result = runMARKStar(make2RL0(), epsilon);

	}

	@Test
	public void test1GUA11() {

		double epsilon = 0.999;
		Stopwatch runtime = new Stopwatch().start();
		Result result = runMARKStar(make1GUA11(), epsilon);
		runtime.stop();

		for (int index = 0; index <6; index++){
			printSequence(result, index);
		}
		// check the results (values collected with e = 0.1 and 64 digits precision)
		assertSequence(result,   0, "HIE VAL", 1.194026e+42, 2.932628e+07, 1.121625e+66, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [7.467257 , 7.467257] (log10)                    complex [66.049848,66.051195] (log10)                    K* = 16.505577 in [16.505576,16.506925] (log10)
		assertSequence(result,   1, "HIE HID", 1.194026e+42, 5.738568e+07, 3.346334e+66, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [7.758803 , 7.758803] (log10)                    complex [66.524569,66.543073] (log10)                    K* = 16.688752 in [16.688752,16.707256] (log10)
		assertSequence(result,   2, "HIE HIE", 1.194026e+42, 6.339230e+06, 5.544100e+65, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [6.802036 , 6.802036] (log10)                    complex [65.743831,65.769366] (log10)                    K* = 16.864781 in [16.864780,16.890316] (log10)
		assertSequence(result,   3, "HIE LYS", 1.194026e+42, 6.624443e+04, 3.315130e+63, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [4.821149 , 4.826752] (log10)                    complex [63.520501,63.563549] (log10)                    K* = 16.622337 in [16.616735,16.665386] (log10)
		assertSequence(result,   4, "HIE ARG", 1.194026e+42, 1.196619e+05, 5.375633e+64, epsilon); // protein [42.077014,42.077014] (log10)                    ligand [5.077956 , 5.087238] (log10)                    complex [64.730430,64.774106] (log10)                    K* = 17.575460 in [17.566178,17.619136] (log10)
		assertSequence(result,   5, "HID VAL", 9.813429e+41, 2.932628e+07, 2.680104e+66, epsilon); // protein [41.991821,41.992159] (log10)                    ligand [7.467257 , 7.467257] (log10)                    complex [66.428152,66.446408] (log10)                    K* = 16.969074 in [16.968735,16.987330] (log10)
		System.out.println("Total time: "+runtime.getTime(2));
	}
	public static void printSequence(Result result, int sequenceIndex){
		MARKStar.ScoredSequence scoredSequence =result.scores.get(sequenceIndex);
		String out = "Printing sequence "+sequenceIndex+": "+scoredSequence.sequence.toString(Sequence.Renderer.ResType)+"\n"+
				"Protein LB: "+String.format("%6.3e",scoredSequence.score.protein.values.qstar)+
				" Protein UB: "+String.format("%6.3e",scoredSequence.score.protein.values.pstar)+"\n"+
				"Ligand LB: "+String.format("%6.3e",scoredSequence.score.ligand.values.qstar)+
				" Ligand UB: "+String.format("%6.3e",scoredSequence.score.ligand.values.pstar)+"\n"+
				"Complex LB: "+String.format("%6.3e",scoredSequence.score.complex.values.qstar)+
				" Complex UB: "+String.format("%6.3e",scoredSequence.score.complex.values.pstar)+"\n"+
				"KStar Score: "+String.format("%6.3e",scoredSequence.score.complex.values.pstar.divide(scoredSequence.score.ligand.values.qstar.multiply(scoredSequence.score.protein.values.qstar), RoundingMode.HALF_UP));
		System.out.println(out);
	}

	public static void assertSequence(Result result, int sequenceIndex, String sequence, Double proteinQStar, Double ligandQStar, Double complexQStar, double epsilon) {

		MARKStar.ScoredSequence scoredSequence = result.scores.get(sequenceIndex);

		// check the sequence
		assertThat(scoredSequence.sequence.toString(Sequence.Renderer.ResType), is(sequence));

		// check q* values and epsilon
		assertResult(scoredSequence.score.protein, proteinQStar, epsilon);
		assertResult(scoredSequence.score.ligand, ligandQStar, epsilon);
		assertResult(scoredSequence.score.complex, complexQStar, epsilon);
	}

	public static void assertResult(PartitionFunction.Result result, Double qstar, double epsilon) {
		if (qstar != null) {
			double comparison = qstar*(1.0 - epsilon);
			if(comparison < 1e-5)
				comparison = 0;
			assertThat(result.status, is(PartitionFunction.Status.Estimated));
			assertThat(result.values.qstar.doubleValue(), greaterThanOrEqualTo(comparison));
			assertThat(result.values.getEffectiveEpsilon(), lessThanOrEqualTo(epsilon));
		} else {
			assertThat(result.status, is(not(PartitionFunction.Status.Estimated)));
		}
	}

	public static ConfSpaces make1A0R() {

		ConfSpaces confSpaces = new ConfSpaces();

		// configure the forcefield
		confSpaces.ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.read(FileTools.readResource("/1A0R_prepped.pdb"));

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
				.addMoleculeForWildTypeRotamers(mol)
				.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("039", "0339")
				.build();
		protein.flexibility.get("0313").setLibraryRotamers(Strand.WildType, "ASN", "SER", "THR", "GLN", "HID", "ALA", "VAL", "ILE", "LEU", "GLY", "ALA", "VAL", "GLY").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("0311").setLibraryRotamers(Strand.WildType, "HID", "HIE", "ASP", "GLU", "SER", "THR", "ASN", "GLN", "ALA", "VAL", "ILE", "LEU", "GLY").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("0332").setLibraryRotamers(Strand.WildType, "TRP", "ALA", "VAL", "ILE", "LEU", "PHE", "TYR", "MET", "SER", "THR", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("0520", "0729")
				.build();
		ligand.flexibility.get("0605").setLibraryRotamers(Strand.WildType, "LEU", "ILE", "ALA", "VAL", "PHE", "TYR", "MET", "GLU", "ASP", "HID", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("0696").setLibraryRotamers(Strand.WildType, "GLU", "ASP", "PHE", "TYR", "ALA", "VAL", "ILE", "LEU", "HIE", "HID", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("0697").setLibraryRotamers(Strand.WildType, "LEU", "ILE", "ALA", "VAL", "PHE", "TYR", "MET", "GLU", "ASP", "HID", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("0698").setLibraryRotamers(Strand.WildType, "LEU", "ILE", "ALA", "VAL", "PHE", "TYR", "MET", "GLU", "ASP", "HID", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("0601").setLibraryRotamers(Strand.WildType, "MET", "ILE", "ALA", "VAL", "LEU", "PHE", "TYR", "GLU", "ASP", "HID", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("0729").setLibraryRotamers(Strand.WildType, "GLU", "ASP", "PHE", "TYR", "ALA", "VAL", "ILE", "LEU", "HIE", "HID", "ASN", "GLN", "GLY").addWildTypeRotamers().setContinuous();

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
	private static List<MARKStar.ScoredSequence> debugFatal() {
		ConfSpaces confSpaces = make1A0R();
		Parallelism parallelism = Parallelism.makeCpu(NUM_CPUs);

		// Define the minimizing energy calculator
		EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
				.setParallelism(parallelism)
				.setIsMinimizing(true)
				.build();
		// Define the rigid energy calculator
		EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
				.setIsMinimizing(false)
				.build();
		// how should we define energies of conformations?
		MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
			return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
					.setEnergyPartition(EnergyPartition.AllOnPairs)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, minimizingEcalc)
							.build()
							.calcReferenceEnergies()
					)
					.build();
		};

		MARKStar.Settings settings = new MARKStar.Settings.Builder()
				.setEpsilon(0.1)
				.setEnergyMatrixCachePattern("*testmat.emat")
				.setShowPfuncProgress(true)
				.setParallelism(parallelism)
				.setReduceMinimizations(REUDCE_MINIMIZATIONS)
				.build();
		MARKStar run = new MARKStar(confSpaces.protein, confSpaces.ligand,
				confSpaces.complex, rigidEcalc, minimizingEcalc, confEcalcFactory, settings);
		return run.run();
	}

	@Test
	public void testFatalError() {
		List<MARKStar.ScoredSequence> markStarSeqs = debugFatal();
		for(MARKStar.ScoredSequence seq: markStarSeqs)
			printMARKStarComputationStats(seq);
	}
}


