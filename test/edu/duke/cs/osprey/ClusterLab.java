package edu.duke.cs.osprey;

import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import org.slf4j.bridge.SLF4JBridgeHandler;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static edu.duke.cs.osprey.tools.Log.log;


public class ClusterLab {

	public static void main(String[] args)
	throws Exception {

		// configure logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();

		// if this is a fork, jump to the fork code
		String idStr = System.getProperty("fork.id");
		if (idStr != null) {

			// get the fork id and size
			int id = Integer.parseInt(idStr);
			int size = Integer.parseInt(System.getProperty("fork.size"));

			run(id, size);
			return;
		}

		// not a fork, so do the forking

		// fork into multiple processes
		// (please don't fork bomb. please, please, please...)
		log("MAIN: forking ...");
		int numForks = 2;
		List<Fork> forks = IntStream.range(0, numForks)
			.mapToObj(id -> new Fork(id, numForks))
			.collect(Collectors.toList());
		log("MAIN: forked!");

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

	private static void run(int id, int size) {

		Parallelism parallelism = Parallelism.makeCpu(4);
		Cluster cluster = new Cluster("Osprey", id, size, parallelism, true);

		// set up a toy design
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();

		if (id > 0) {

			// run as a pure member node
			try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams).build()) {

				Cluster.Member member = cluster.new Member();
				member.putTaskContext(EnergyTask.class, new EnergyTask.Context(confSpaces.complex, ecalc));
				member.run();
			}

		} else {

			// run as a client/member node
			try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
				.setCluster(cluster)
				.build()) {

				Cluster.Client tasks = (Cluster.Client)ecalc.tasks;
				tasks.putContext(EnergyTask.class, new EnergyTask.Context(confSpaces.complex, ecalc));

				int[][] confs = {
					{ 0, 0, 0, 0, 0, 0, 0, 0 },
					{ 1, 0, 0, 0, 0, 0, 0, 0 },
					{ 0, 1, 0, 0, 0, 0, 0, 0 },
					{ 0, 0, 1, 0, 0, 0, 0, 0 },
					{ 0, 0, 0, 1, 0, 0, 0, 0 },
					{ 0, 0, 0, 0, 1, 0, 0, 0 },
					{ 0, 0, 0, 0, 0, 1, 0, 0 },
					{ 0, 0, 0, 0, 0, 0, 1, 0 },
					{ 0, 0, 0, 0, 0, 0, 0, 1 }
				};
				for (int i=0; i<100; i++) {
					final int fi = i;
					tasks.submit(
						new EnergyTask(i, confs[i % confs.length]),
						energy -> log("CLIENT: conf %d energy = %f", fi, energy)
					);
				}
				tasks.waitForFinish();
			}
		}

		/* TODO: NEXTTIME: try to get K* working
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

			SimpleConfSpace confSpace = (SimpleConfSpace)info.confSpace;

			// how should we define energies of conformations?
			info.confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
					.build()
					.calcReferenceEnergies()
				)
				.build();

			// calc energy matrix
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
				.build()
				.calcEnergyMatrix();

			// how should confs be ordered and searched?
			info.confSearchFactory = (rcs) ->
				new ConfAStarTree.Builder(emat, rcs)
					.setTraditional()
					.build();
		}

		// run K*
		TestKStar.Result result = new TestKStar.Result();
		result.kstar = kstar;
		result.scores = kstar.run();
		*/

		if (id > 0) {
			log("MEMBER %d: finished", id);
		} else {
			log("CLIENT: finished");
		}
	}
}


class EnergyTask extends Cluster.Task<Double,EnergyTask.Context> {

	public static class Context {

		final SimpleConfSpace confSpace;
		final EnergyCalculator ecalc;

		public Context(SimpleConfSpace confSpace, EnergyCalculator ecalc) {
			this.confSpace = confSpace;
			this.ecalc = ecalc;
		}
	}

	final int i;
	final int[] conf;

	EnergyTask(int i, int[] conf) {
		this.i = i;
		this.conf = conf;
	}

	@Override
	public Double run(Context ctx) {

		RCTuple frag = new RCTuple(conf);
		ParametricMolecule pmol = ctx.confSpace.makeMolecule(conf);
		ResidueInteractions inters = ResInterGen.of(ctx.confSpace)
			.addIntras(frag)
			.addInters(frag)
			.addShell(frag)
			.make();

		return ctx.ecalc.calcEnergy(pmol, inters).energy;
	}
}
