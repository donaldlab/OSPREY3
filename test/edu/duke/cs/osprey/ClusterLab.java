package edu.duke.cs.osprey;

import edu.duke.cs.osprey.parallelism.Cluster;
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

	private static final String ClusterName = "ClusterLab";

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

			if (id == 0) {
				client(id, size);
			} else {
				member(id, size);
			}
			return;
		}

		// not a fork, so do the forking

		// fork into multiple processes
		// (please don't fork bomb. please, please, please...)
		log("MAIN: forking ...");
		int numForks = 3;
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

	private static void member(int id, int size) {
		new Cluster.Member(ClusterName, id, size).run();
	}

	private static void client(int id, int size) {

		log("CLIENT: start");

		try (Cluster.Client client = new Cluster.Client(ClusterName, id, size)) {

			log("CLIENT: ready");

			// damn, Java is so verbose sometimes...
			class FooTask extends Cluster.Task<String> {

				final int i;

				FooTask(int i) {
					this.i = i;
				}

				@Override
				public String run() {
					try { Thread.sleep(2000); } catch (Throwable t) {} // TEMP
					return "Foo" + i;
				}
			}

			try (TaskExecutor tasks = client.tasks()) {

				log("CLIENT: sending tasks to %d nodes ...", tasks.getParallelism());
				for (int i=0; i<10; i++) {
					tasks.submit(
						new FooTask(i),
						(result) -> {
							log("CLIENT: result = %s", result);
						}
					);
					log("CLIENT: submitted %d", i);
				}
				log("CLIENT: all submitted, waiting ...");
				tasks.waitForFinish();
			}
		}

		log("CLIENT: finished");
	}
}
