package edu.duke.cs.osprey.parallelism;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static edu.duke.cs.osprey.tools.Log.log;


public class ForkCluster {

	public static void run(int numNodes, boolean clientIsMember, Class<?> mainClass, Consumer<Cluster> task) {

		final String clusterName = "ForkCluster";
		final String jobId = "fork";

		// if this is a fork, jump to the fork code
		String idStr = System.getProperty("fork.id");
		if (idStr != null) {

			// get the fork id and size
			int id = Integer.parseInt(idStr);

			task.accept(new Cluster(clusterName, jobId, id, numNodes, clientIsMember));
			return;
		}

		// not a fork, so do the forking

		// fork into multiple processes
		// (please don't fork bomb. please, please, please...)
		log("MAIN: forking ...");
		List<Process> forks = IntStream.range(1, numNodes)
			.mapToObj(id -> {

				// start the JVM process
				ProcessBuilder pb = new ProcessBuilder(
					Paths.get(System.getProperty("java.home")).resolve("bin").resolve("java").toString(),
					"-cp", System.getProperty("java.class.path"),
					"-Dfork.id=" + id,
					mainClass.getCanonicalName()
				);
				pb.directory(new File(System.getProperty("user.dir")));
				pb.redirectOutput(ProcessBuilder.Redirect.INHERIT);
				pb.redirectError(ProcessBuilder.Redirect.INHERIT);

				try {
					return pb.start();
				} catch (IOException ex) {
					throw new RuntimeException("can't start JVM process", ex);
				}
			})
			.collect(Collectors.toList());
		log("MAIN: forked!");

		// run the client here, so we can cancel it from the IDE
		task.accept(new Cluster(clusterName, jobId, 0, numNodes, clientIsMember));

		// wait for the forks to finish
		for (Process fork : forks) {
			try {
				fork.waitFor();
			} catch (InterruptedException ex) {
				ex.printStackTrace(System.err);
			}
		}
		log("MAIN: all forks complete!");
	}
}
