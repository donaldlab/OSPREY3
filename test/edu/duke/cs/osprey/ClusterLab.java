package edu.duke.cs.osprey;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.slf4j.bridge.SLF4JBridgeHandler;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static edu.duke.cs.osprey.tools.Log.log;


public class ClusterLab {

	enum SourceConfSpace {
		Test2RL0,
		TestCOVIDSmall,
		TestCOVIDActual,

	}

	public static void main(String[] args)
	throws Exception {

		// configure logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();

		//forkCluster();
		//multiProcessCluster(args);
		slurmCluster();
	}

	private static void forkCluster()
	throws Exception {

		final String clusterName = "ForkCluster";
		final String jobId = "fork";
		final int numNodes = 2;
		final boolean clientIsMember = false;

		// if this is a fork, jump to the fork code
		String idStr = System.getProperty("fork.id");
		if (idStr != null) {

			// get the fork id and size
			int id = Integer.parseInt(idStr);
			int size = Integer.parseInt(System.getProperty("fork.size"));

			run(new Cluster(clusterName, jobId, id, size, clientIsMember));
			return;
		}

		// not a fork, so do the forking

		// fork into multiple processes
		// (please don't fork bomb. please, please, please...)
		log("MAIN: forking ...");
		List<Fork> forks = IntStream.range(1, numNodes)
			.mapToObj(id -> new Fork(id, numNodes))
			.collect(Collectors.toList());
		log("MAIN: forked!");

		// run the client here, so we can cancel it from the IDE
		run(new Cluster(clusterName, jobId, 0, numNodes, clientIsMember));

		// wait for the forks to finish
		for (Fork fork : forks) {
			fork.process.waitFor();
		}
		log("MAIN: all forks complete!");
	}

	private static class Fork {

		final int id;
		final Process process;

		Fork(int id, int size) {

			this.id = id;

			// start the JVM process
			ProcessBuilder pb = new ProcessBuilder(
				Paths.get(System.getProperty("java.home")).resolve("bin").resolve("java").toString(),
				"-cp", System.getProperty("java.class.path"),
				"-Dfork.id=" + id,
				"-Dfork.size=" + size,
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
		int numNodes = Integer.parseInt(args[2]);
		boolean clientIsMember = false;

		run(new Cluster("MultiProcessCluster", jobId, nodeId, numNodes, clientIsMember));
	}

	private static void slurmCluster() {
		run(Cluster.fromSLURM(true));
	}

	private static void run(Cluster cluster) {

		Stopwatch stopwatch = new Stopwatch().start();

		Parallelism parallelism = Parallelism.makeCpu(4);

		// set up a toy design
		TestKStar.ConfSpaces confSpaces = null;
		SourceConfSpace source = SourceConfSpace.Test2RL0;
		switch (source) {
			case Test2RL0:
			// set up a toy design
			confSpaces = TestKStar.make2RL0();
			break;
			case TestCOVIDSmall:
			confSpaces = TestCOVID.makeCOVIDSmall();
			break;
			case TestCOVIDActual:
			confSpaces = TestCOVID.makeMakeCOVIDActual();
		}

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces.asList(), confSpaces.ffparams)
			.setCluster(cluster)
			.setParallelism(parallelism)
			.build()) {

			// run K*
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

			Molecule mol = PDBIO.readResource("test-resources/ACE2_1r4l_peptide5vay_complex_model_leidosprep_declashed4.pdb");

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

		public static TestKStar.ConfSpaces makeMakeCOVIDActual() {
			TestKStar.ConfSpaces confSpaces = new TestKStar.ConfSpaces();

			// configure the forcefield
			confSpaces.ffparams = new ForcefieldParams();

			Molecule mol = PDBIO.readResource("test-resources/ACE2_1r4l_peptide5vay_complex_model_leidosprep_declashed4.pdb");

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
	}
}

