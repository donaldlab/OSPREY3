package edu.duke.cs.osprey;

import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.MALEEC;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.Slurm;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.slf4j.bridge.SLF4JBridgeHandler;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static edu.duke.cs.osprey.tools.Log.log;


public class ClusterLab {

	enum SourceConfSpace {
		Test2RL0,
		TestCOVIDSmall,
		TestCOVIDActual
	}

	public static void main(String[] args)
	throws Exception {

		// configure logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();

		forkCluster();
		//multiProcessCluster(args);
		//slurmCluster();
	}

	private static void forkCluster()
	throws Exception {

		final String clusterName = "ForkCluster";
		final String jobId = "fork";
		final int numNodes = 1;
		List<String> nodes = IntStream.range(0, numNodes)
			.mapToObj(i -> "localhost")
			.collect(Collectors.toList());
		final boolean clientIsMember = true;

		// if this is a fork, jump to the fork code
		String idStr = System.getProperty("fork.id");
		if (idStr != null) {

			// get the fork id and size
			int id = Integer.parseInt(idStr);

			run(new Cluster(clusterName, jobId, id, nodes, clientIsMember));
			return;
		}

		// not a fork, so do the forking

		// fork into multiple processes
		// (please don't fork bomb. please, please, please...)
		log("MAIN: forking ...");
		List<Fork> forks = IntStream.range(1, numNodes)
			.mapToObj(Fork::new)
			.collect(Collectors.toList());
		log("MAIN: forked!");

		// run the client here, so we can cancel it from the IDE
		run(new Cluster(clusterName, jobId, 0, nodes, clientIsMember));

		// wait for the forks to finish
		for (Fork fork : forks) {
			fork.process.waitFor();
		}
		log("MAIN: all forks complete!");
	}

	private static class Fork {

		final int id;
		final Process process;

		Fork(int id) {

			this.id = id;

			// start the JVM process
			ProcessBuilder pb = new ProcessBuilder(
				Paths.get(System.getProperty("java.home")).resolve("bin").resolve("java").toString(),
				"-cp", System.getProperty("java.class.path"),
				"-Dfork.id=" + id,
				ClusterLab.class.getCanonicalName()
			);
			pb.directory(new File(System.getProperty("user.dir")));
			pb.redirectOutput(ProcessBuilder.Redirect.INHERIT);
			pb.redirectError(ProcessBuilder.Redirect.INHERIT);
			try {
				process = pb.start();
			} catch (IOException ex) {
				throw new RuntimeException("can't start JVM process", ex);
			}
		}
	}

	private static void multiProcessCluster(String[] args) {

		// parse the args to get nodeId and cluster size
		String jobId = args[0];
		int nodeId = Integer.parseInt(args[1]);
		boolean clientIsMember = false;
		List<String> nodes = Arrays.asList("localhost");

		run(new Cluster("MultiProcessCluster", jobId, nodeId, nodes, clientIsMember));
	}

	private static void slurmCluster() {
		run(Slurm.makeCluster(true));
	}

	private static void run(Cluster cluster) {

		Stopwatch stopwatch = new Stopwatch().start();

		int threads;
		try {
			// try to read the threads from SLURM, if possible
			threads = Integer.parseInt(System.getenv("SLURM_CPUS_PER_TASK"));
		} catch (Throwable t) {
			threads = 4;
		}
		Parallelism parallelism = Parallelism.makeCpu(threads);

		if (cluster.nodeId > 0) {
			log("MEMBER %d: started with %d threads", cluster.nodeId, threads);
		} else {
			log("CLIENT: started with %d threads", threads);
		}

		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();
		//TestKStar.ConfSpaces confSpaces = TestCOVID.makeCOVIDSmall();
		//TestKStar.ConfSpaces confSpaces = TestCOVID.makeMakeCOVIDComplexMedium();
		//TestKStar.ConfSpaces confSpaces = TestCOVID.makeMakeCOVIDActual();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces.asList(), confSpaces.ffparams)
			.setParallelism(parallelism)
			.build()) {

			/* run K*
				or not. MALEEC to the rescue!!
			KStarScoreWriter.Formatter testFormatter = (KStarScoreWriter.ScoreInfo info) -> {

				Function<PartitionFunction.Result,String> formatPfunc = (pfuncResult) -> {
					if (pfuncResult.status == PartitionFunction.Status.Estimated) {
						return String.format("%12e", pfuncResult.values.qstar.doubleValue());
					}
					return "null";
				};

				return String.format("assertSequence(result, %3d, \"%s\", %-12s, %-12s, %-12s, epsilon); // protein %s ligand %s complex %s K* = %s",
					info.sequenceNumber,
					info.sequence.toString(Sequence.Renderer.ResType),
					formatPfunc.apply(info.kstarScore.protein),
					formatPfunc.apply(info.kstarScore.ligand),
					formatPfunc.apply(info.kstarScore.complex),
					info.kstarScore.protein.toString(),
					info.kstarScore.ligand.toString(),
					info.kstarScore.complex.toString(),
					info.kstarScore.toString()
				);
			};

			// configure K*
			KStar.Settings settings = new KStar.Settings.Builder()
				.setEpsilon(0.95)
				.setStabilityThreshold(null)
				.addScoreConsoleWriter(testFormatter)
				.setMaxSimultaneousMutations(1)
				.build();
			KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, settings);
			for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

				// turn off the default confdb for tests
				info.confDBFile = null;

				SimpleConfSpace confSpace = (SimpleConfSpace)info.confSpace;

				// compute reference energies locally on each node
				SimpleReferenceEnergies eref;
				try (EnergyCalculator localEcalc = ecalc.local()) {
					eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, localEcalc)
						.build()
						.calcReferenceEnergies();
				}

				// how should we define energies of conformations?
				info.confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
					.setReferenceEnergies(eref)
					.build();

				// calc the energy matrix
				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
					.build()
					.calcEnergyMatrix();

				info.pfuncFactory = rcs -> new GradientDescentPfunc(
					info.confEcalc,
					// HACKHACK: don't make A* trees on member nodes, we don't have the energy matrix
					cluster.nodeId > 0 ? null : new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build(),
					cluster.nodeId > 0 ? null : new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build(),
					rcs.getNumConformations()
				);
			}

			// run K*
			kstar.run(ecalc.tasks);
			*/

			// we want the complex conf space, right?
			SimpleConfSpace confSpace = confSpaces.complex;

			// pick your favorite sequence
			Sequence seq = confSpace.seqSpace.makeWildTypeSequence()
				.set("G649", "TYR");

			if (!seq.isFullyAssigned()) {
				throw new Error();
			}

			// compute reference energies locally on each node
			SimpleReferenceEnergies eref;
			try (EnergyCalculator localEcalc = ecalc.local()) {
				eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, localEcalc)
					.build()
					.calcReferenceEnergies();
			}

			// how should we define energies of conformations?
			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				.setReferenceEnergies(eref)
				.setEnergyPartition(EnergyPartition.AllOnPairs) // get tigher lower bounds on energies!
				.build();

			// calc the energy matrix
			// PROTIP: If your conformation space has only one sequence in it,
			// emat calculation will go waaaay faster
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(new File("maleec.emat"))
				//.setTripleCorrectionThreshold(1.0)
				.build()
				.calcEnergyMatrix();

			// massively awesome low-energy ensemble calculator
			MALEEC maleec = new MALEEC(ecalc.tasks, confEcalc, emat);
			maleec.doEeet(seq, 50, 60*10 /* 10 minutes */, "");
		}

		if (cluster.nodeId > 0) {
			log("MEMBER %d: finished in %s", cluster.nodeId, stopwatch.stop().getTime(2));
		} else {
			log("CLIENT: finished in %s", stopwatch.stop().getTime(2));
		}
	}

	private static class TestCOVID {

		public static TestKStar.ConfSpaces makeCOVIDSmall() {
			TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();

			// configure the forcefield
			confSpaces.ffparams = new ForcefieldParams();

			Molecule mol = PDBIO.readResource("/ACE2_1r4l_peptide5vay_complex_model_leidosprep_declashed4.pdb");

			// make sure all strands share the same template library
			ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

			// define the protein strand
			Strand protein = new Strand.Builder(mol)
					.setTemplateLibrary(templateLib)
					.setResidues("B123", "B139")
					.build();
			protein.flexibility.get("B125").setLibraryRotamers(Strand.WildType, "ARG", "TRP", "PHE", "TYR", "LYS", "HIP", "HIE", "HID", "MET", "ASP", "ASN", "GLU", "GLN").addWildTypeRotamers().setContinuous();
			protein.flexibility.get("B129").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

			// define the ligand strand
			Strand ligand = new Strand.Builder(mol)
					.setTemplateLibrary(templateLib)
					.setResidues("A19", "A615")
					.build();
			ligand.flexibility.get("A329").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A330").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

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

		public static TestKStar.ConfSpaces makeMakeCOVIDComplexMedium() {
			TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();

			// configure the forcefield
			confSpaces.ffparams = new ForcefieldParams();

			Molecule mol = PDBIO.readResource("/ACE2_1r4l_peptide5vay_complex_model_leidosprep_declashed4.pdb");

			// make sure all strands share the same template library
			ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

			// define the protein strand
			Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("B123", "B139")
				.build();
			// mutable
			protein.flexibility.get("B123").setLibraryRotamers("ARG").setContinuous(); // wt=ASN
			protein.flexibility.get("B125").setLibraryRotamers("HIP").setContinuous(); // wt=THR
			protein.flexibility.get("B126").setLibraryRotamers("HIP").setContinuous(); // wt=ALA
			protein.flexibility.get("B129").setLibraryRotamers("ARG").setContinuous(); // wt=ARG *same*
			protein.flexibility.get("B130").setLibraryRotamers("ARG").setContinuous(); // wt=PHE
			protein.flexibility.get("B133").setLibraryRotamers("THR").setContinuous(); // wt=ALA
			protein.flexibility.get("B134").setLibraryRotamers("HIP").setContinuous(); // wt=VAL
			protein.flexibility.get("B136").setLibraryRotamers("ARG").setContinuous(); // wt=GLU
			protein.flexibility.get("B138").setLibraryRotamers("ARG").setContinuous(); // wt=PHE
			// flexible
			//protein.flexibility.get("B124").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			//protein.flexibility.get("B139").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

			// define the ligand strand
			Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A19", "A615")
				.build();
			ligand.flexibility.get("A34").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A35").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A37").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A38").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			//ligand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			//ligand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			//ligand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			//ligand.flexibility.get("A49").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			//ligand.flexibility.get("A329").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			//ligand.flexibility.get("A330").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			//ligand.flexibility.get("A352").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			//ligand.flexibility.get("A353").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			//ligand.flexibility.get("A354").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			//ligand.flexibility.get("A355").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			//ligand.flexibility.get("A357").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

			// make the conf spaces ("complex" SimpleConfSpace, har har!)
			confSpaces.complex = new SimpleConfSpace.Builder()
				.addStrands(protein, ligand)
				.build();

			// the others can't be left null, so just make it complexes all the way down
			confSpaces.protein = confSpaces.complex;
			confSpaces.ligand = confSpaces.complex;

			return confSpaces;
		}

		public static TestKStar.ConfSpaces makeMakeCOVIDActual() {
			TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();

			// configure the forcefield
			confSpaces.ffparams = new ForcefieldParams();

			Molecule mol = PDBIO.readResource("/ACE2_1r4l_peptide5vay_complex_model_leidosprep_declashed4.pdb");

			// make sure all strands share the same template library
			ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

			// define the protein strand
			Strand protein = new Strand.Builder(mol)
					.setTemplateLibrary(templateLib)
					.setResidues("B123", "B139")
					.build();
			// mutable
			protein.flexibility.get("B123").setLibraryRotamers("ARG").setContinuous();
			protein.flexibility.get("B125").setLibraryRotamers("HIP").setContinuous();
			protein.flexibility.get("B126").setLibraryRotamers("HIP").setContinuous();
			protein.flexibility.get("B129").setLibraryRotamers("ARG").setContinuous();
			protein.flexibility.get("B130").setLibraryRotamers("ARG").setContinuous();
			protein.flexibility.get("B133").setLibraryRotamers("THR").setContinuous();
			protein.flexibility.get("B134").setLibraryRotamers("HIP").setContinuous();
			protein.flexibility.get("B136").setLibraryRotamers("ARG").setContinuous();
			protein.flexibility.get("B138").setLibraryRotamers("ARG").setContinuous();
			// flexible
			protein.flexibility.get("B124").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			protein.flexibility.get("B139").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

			// define the ligand strand
			Strand ligand = new Strand.Builder(mol)
					.setTemplateLibrary(templateLib)
					.setResidues("A19", "A615")
					.build();
			ligand.flexibility.get("A34").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A35").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A37").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A38").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A49").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A329").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A330").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A352").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A353").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A354").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A355").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A357").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

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

		public static TestKStar.ConfSpaces makeMakeCOVIDSmallerMonster() {
			TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();

			// configure the forcefield
			confSpaces.ffparams = new ForcefieldParams();

			Molecule mol = PDBIO.readResource("/ACE2_1r42_peptide3qf7_complex_model_leidosprep_declashed_trimmed3.pdb");

			// make sure all strands share the same template library
			ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

			// define the protein strand
			Strand protein = new Strand.Builder(mol)
					.setTemplateLibrary(templateLib)
					.setResidues("B166", "B190")
					.build();
			// mutable
			protein.flexibility.get("B166").setLibraryRotamers("SER").setContinuous();
			protein.flexibility.get("B169").setLibraryRotamers("ARG").setContinuous();
			protein.flexibility.get("B170").setLibraryRotamers("GLN").setContinuous();
			protein.flexibility.get("B173").setLibraryRotamers("HIP").setContinuous();
			protein.flexibility.get("B176").setLibraryRotamers("ARG").setContinuous();
			protein.flexibility.get("B177").setLibraryRotamers("ASN").setContinuous();
			protein.flexibility.get("B180").setLibraryRotamers("LYS").setContinuous();
			protein.flexibility.get("B181").setLibraryRotamers("TYR").setContinuous();
			protein.flexibility.get("B184").setLibraryRotamers("GLN").setContinuous();
			protein.flexibility.get("B185").setLibraryRotamers("GLU").setContinuous();
			protein.flexibility.get("B187").setLibraryRotamers("ARG").setContinuous();
			protein.flexibility.get("B188").setLibraryRotamers("ARG").setContinuous();
			// flexible
			protein.flexibility.get("B167").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

			// define the ligand strand
			Strand ligand = new Strand.Builder(mol)
					.setTemplateLibrary(templateLib)
					.setResidues("A19", "A615")
					.build();
			setMonsterFlexFirstForm(ligand);

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

		public static TestKStar.ConfSpaces makeMakeCOVIDBossMonster() {
			TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();

			// configure the forcefield
			confSpaces.ffparams = new ForcefieldParams();

			Molecule mol = PDBIO.readResource("/ACE2_1r42_peptide3qf7_complex_model_leidosprep_declashed_trimmed3.pdb");

			// make sure all strands share the same template library
			ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld).build();

			// define the protein strand
			Strand protein = new Strand.Builder(mol)
					.setTemplateLibrary(templateLib)
					.setResidues("B166", "B190")
					.build();
			setMonterSeq1(protein);
			// flexible
			protein.flexibility.get("B168").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			protein.flexibility.get("B171").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			protein.flexibility.get("B172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			protein.flexibility.get("B174").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			protein.flexibility.get("B178").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			protein.flexibility.get("B183").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();


			// define the ligand strand
			Strand ligand = new Strand.Builder(mol)
					.setTemplateLibrary(templateLib)
					.setResidues("A19", "A615")
					.build();
			setMonsterFlexFirstForm(ligand);

			// Second form revealed!
			ligand.flexibility.get("A53").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A57").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A58").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A64").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A68").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A334").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

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

		private static void setMonsterFlexFirstForm(Strand ligand) {
			ligand.flexibility.get("A35").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A38").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A39").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A46").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A48").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A49").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A52").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A61").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A329").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A330").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A331").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A332").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A340").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A353").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A355").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
			ligand.flexibility.get("A357").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		}

		private static void setMonterSeq1(Strand protein) {
			// mutable
			protein.flexibility.get("B166").setLibraryRotamers("SER").setContinuous();
			protein.flexibility.get("B169").setLibraryRotamers("ARG").setContinuous();
			protein.flexibility.get("B170").setLibraryRotamers("GLN").setContinuous();
			protein.flexibility.get("B173").setLibraryRotamers("HIP").setContinuous();
			protein.flexibility.get("B176").setLibraryRotamers("ARG").setContinuous();
			protein.flexibility.get("B177").setLibraryRotamers("ASN").setContinuous();
			protein.flexibility.get("B180").setLibraryRotamers("LYS").setContinuous();
			protein.flexibility.get("B181").setLibraryRotamers("TYR").setContinuous();
			protein.flexibility.get("B184").setLibraryRotamers("GLN").setContinuous();
			protein.flexibility.get("B185").setLibraryRotamers("GLU").setContinuous();
			protein.flexibility.get("B187").setLibraryRotamers("ARG").setContinuous();
			protein.flexibility.get("B188").setLibraryRotamers("ARG").setContinuous();
		}
	}
}

