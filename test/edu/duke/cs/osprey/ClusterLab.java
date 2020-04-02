package edu.duke.cs.osprey;

import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import org.slf4j.bridge.SLF4JBridgeHandler;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;

import static edu.duke.cs.osprey.tools.Log.log;


public class ClusterLab {

	private static final String ClusterName = "ClusterLab";

	public static void main(String[] args)
	throws Exception {

		// configure logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();

		// if this is a fork, jump to the fork code
		String idStr = System.getProperty("fork.id");
		if (idStr != null) {
			int id = Integer.parseInt(idStr);
			if (id == 0) {
				client(id);
			} else {
				member(id);
			}
			return;
		}

		// not a fork, so do the forking

		// fork into multiple processes
		// (please don't fork bomb. please, please, please...)
		log("MAIN: forking ...");
		Fork[] forks = {
			new Fork(0),
			new Fork(1),
			new Fork(2)
		};
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

		Fork(int id)
		throws IOException {

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
			process = pb.start();
		}
	}

	private static void member(int id) {
		new Cluster.Member(ClusterName, id).run();
	}

	private static void client(int id) {

		log("CLIENT: start");

		try (Cluster.Client client = new Cluster.Client(ClusterName, id)) {

			log("CLIENT: ready");

			// damn, Java is so verbose sometimes...
			class FooTask extends Cluster.Task<String> {

				final int i;

				FooTask(int i) {
					this.i = i;
				}

				@Override
				public String run() {
					return "Foo" + i;
				}
			}

			TaskExecutor tasks = client.randomTasks();
			log("CLIENT: sending tasks to %d nodes ...", tasks.getParallelism());
			for (int i=0; i<3; i++) {
				tasks.submit(
					new FooTask(i),
					(result) -> {
						log("CLIENT: result = %s", result);
					}
				);
			}
			log("CLIENT: submitted, waiting ...");
			tasks.waitForFinish();
		}

		log("CLIENT: finished");
	}
}
