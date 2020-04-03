package edu.duke.cs.osprey;

import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
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

		// TEMP
		Parallelism parallelism = new Parallelism(
			4, 0, 0,
			new Parallelism.ClusterInfo("Osprey", id, size, false)
		);
		//Parallelism parallelism = Parallelism.makeCpu(4);

		// set up a toy design
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();

		// how should we compute energies of molecules?
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
			.setParallelism(parallelism)
			.build()) {

			EnergyTask.Context ctx = new EnergyTask.Context(confSpaces.complex, ecalc);

			if (id > 0) {

				// run as a pure member node
				Cluster.Member member = new Cluster.Member(parallelism.clusterInfo);
				member.putTaskContext(EnergyTask.class, ctx);
				member.run(parallelism.getParallelism());

			} else {

				// run as a client/member node
				TaskExecutor.WithContext tasks = (TaskExecutor.WithContext)ecalc.tasks;
				tasks.putContext(EnergyTask.class, ctx);

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

				log("CLIENT: K* finished!");
			}
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
