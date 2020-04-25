package edu.duke.cs.osprey.parallelism;


import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class Slurm {

	/** Read environment variables to determine cluster properties */
	public static Cluster makeCluster(boolean clientIsMember) {
		return new Cluster(
			"SLURM",
			getenv("SLURM_JOB_ID", "SLURM_JOBID"),
			Integer.parseInt(getenv("SLURM_PROCID")),
			parseNodes(
				Integer.parseInt(getenv("SLURM_JOB_NUM_NODES", "SLURM_NNODES")),
				getenv("SLURM_JOB_NODELIST", "SLURM_NODELIST"),
				Integer.parseInt(getenv("SLURM_NTASKS", "SLURM_NPROCS")),
				getenv("SLURM_TASKS_PER_NODE")
			),
			clientIsMember
		);
	}

	private static String getenv(String ... vars) {
		return Arrays.stream(vars)
			.map(System::getenv)
			.filter(Objects::nonNull)
			.findFirst()
			.orElseThrow(() -> new IllegalStateException(String.format("missing envvar: one of %s", Arrays.toString(vars))));
	}

	public static List<String> parseNodes(int numNodes, String nodelist, int numTasks, String tasksPerNode) {

		List<String> nodes = parseNodelist(nodelist);
		if (nodes.size() != numNodes) {
			throw new SlurmException("nodelist " + nodes + " doesn't match number of nodes: " + numNodes);
		}

		List<Integer> tasks = parseTasklist(tasksPerNode);
		if (tasks.stream().mapToInt(i -> i).sum() != numTasks) {
			throw new SlurmException("tasklist " + tasks + " doesn't match number of tasks: " + numTasks);
		}

		if (nodes.size() != tasks.size()) {
			throw new SlurmException("Nodes don't match tasks\nnodes: " + nodes + "\ntasks: " + tasks);
		}

		// copy each node by the number of tasks
		return IntStream.range(0, nodes.size())
			.boxed()
			.flatMap(i ->
				IntStream.range(0, tasks.get(i))
					.mapToObj(t -> nodes.get(i))
			)
			.collect(Collectors.toList());
	}

	private static List<String> parseNodelist(String nodelist) {

		// can be eg:
		// linux[33-34,45-46]
		// linux[38,40-41,47]
		// grisman-39
		// gpu-compute7
		// jerry1
		// a,b,c

		// but not:
		// a[4-5],b[4-5]
		// right?

		// check for simple list first
		int open = nodelist.indexOf('[');
		int close = nodelist.indexOf(']');
		if (open < 0 || close < 0 || close <= open) {

			// easy peasy
			return Arrays.stream(nodelist.split(","))
				.collect(Collectors.toList());
		}

		// ok, now deal with the ranges
		String prefix = nodelist.substring(0, open);
		String content = nodelist.substring(open + 1, close);
		return Arrays.stream(content.split(","))
			.flatMap(val -> {

				// val will be eg:
				// 45
				// apple
				// 9-32

				// is this a range?
				int dash = val.indexOf('-');
				if (dash < 0) {

					// nope, easy peasy
					return Stream.of(prefix + val);
				}

				// hard case, handle the range
				String start = val.substring(0, dash);
				String stop = val.substring(dash + 1);
				try {
					int istart = Integer.parseInt(start);
					int istop = Integer.parseInt(stop);
					return IntStream.rangeClosed(istart, istop)
						.mapToObj(i -> prefix + i);
				} catch (NumberFormatException ex) {
					throw new SlurmException("don't know how to parse nodelist segment: " + val);
				}
			})
			.collect(Collectors.toList());
	}

	private static List<Integer> parseTasklist(String tasksPerNode) {

		// can be eg:
		// 5
		// 1(x3)
		// 1,3,1

		return Arrays.stream(tasksPerNode.split(","))
			.flatMap(section -> {

				try {
					// check for the multiplicative format
					int start = section.indexOf('(');
					int end = section.indexOf(')');
					if (start >= 0 && end > start) {
						int numTasks = Integer.parseInt(section.substring(0, start));
						int numTimes = Integer.parseInt(section.substring(start + 2, end));
						return IntStream.range(0, numTimes)
							.mapToObj(i -> numTasks);
					}

					// otherwise, it's just a plain number
					return Stream.of(Integer.parseInt(section));

				} catch (NumberFormatException ex) {
				throw new SlurmException("don't know how to parse tasklist section: " + section);
				}
			})
			.collect(Collectors.toList());
	}

	public static class SlurmException extends RuntimeException {
		public SlurmException(String msg) {
			super(msg + "\n\t"
				+ System.getenv().entrySet().stream()
					.filter(e -> e.getKey().startsWith("SLURM_"))
					.map(e -> String.format("%20s = %s", e.getKey(), e.getValue()))
					.collect(Collectors.joining("\n\t"))
			);
		}
	}
}
