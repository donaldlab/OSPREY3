package edu.duke.cs.osprey.coffee.nodedb;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.BigExp;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;


public class NodePerformance {

	public static final long InitialNs = 100_000L; // start with 100 us per node operation
	public static final int HistorySize = 100;

	private static class TypePerf {

		Deque<Long> nsHistory = new ArrayDeque<>(HistorySize);
		long nsSum;
		long avgNs;

		TypePerf() {
			clear();
		}

		void clear() {

			nsSum = 0L;
			avgNs = InitialNs;

			// init the history with something neutral, so we can ease-in changes over time
			nsHistory.clear();
			for (int i=0; i<HistorySize; i++) {
				nsHistory.add(InitialNs);
				nsSum += InitialNs;
			}
		}

		void update(long ns) {

			nsSum -= nsHistory.removeFirst();
			nsHistory.add(ns);
			nsSum += ns;
			avgNs = nsSum/HistorySize;

			// just in case
			if (nsSum < 0) {
				throw new Error("node performance should not predict negative node processing times, we don't know how to time travel yet");
			}
		}

		public BigExp score(BigExp zSumUpper) {
			var score = new BigExp(zSumUpper);
			score.div(avgNs);
			score.normalize(true);
			return score;
		}
	}

	private static class StatePerf {

		/** indexed by number assigned positions */
		final TypePerf[] types;

		StatePerf(MultiStateConfSpace.State state) {
			types = new TypePerf[state.confSpace.numPos() + 1];
			for (int i=0; i<types.length; i++) {
				types[i] = new TypePerf();
			}
		}

		void clear() {
			for (var type : types) {
				type.clear();
			}
		}

		TypePerf get(int[] conf) {
			return types[Conf.countAssignments(conf)];
		}
	}

	private class PerfLog {

		final File file;
		final List<String> buf = new ArrayList<>();

		long flushNs = -1;
		boolean failed = false;

		PerfLog(File file) {

			this.file = file;

			// init the log file with a header
			try (var out = new FileWriter(file, false)) {
				out.write(String.join("\t", List.of(
					"state",
					"num assignments",
					"avg Ns",
					"actual Ns",
					"zSumUpper",
					"score",
					"reduction",
					"reduction/Ns",
					"zSumLower",
					"zSumUpperFactor",
					"zSumLowerFactor"
				)));
				out.write("\n");
			} catch (IOException ex) {
				if (!failed) {
					ex.printStackTrace(System.err);
					failed = true;
				}
			}
		}

		void add(int statei, int[] conf, long ns, BigExp zSumUpper, BigExp reduction, BigExp score, BigExp zSumLower, BigExp zSumUpperFactor, BigExp zSumLowerFactor) {

			// compute the actual reduction per unit time
			BigExp actual = new BigExp(reduction);
			actual.div(ns);
			actual.normalize(true);

			// get the current estimation of the average time
			long avgNs = states[statei].get(conf).avgNs;

			// buffer the log entry for efficient file IO
			buf.add(String.format("%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
				statei,
				Conf.countAssignments(conf),
				avgNs,
				ns,
				zSumUpper.log(),
				score.log(),
				reduction.log(),
				actual.log(),
				zSumLower.log(),
				zSumUpperFactor.log(),
				zSumLowerFactor.log()
			));

			// flush to disk at most every second
			if (System.nanoTime() - flushNs >= 1_000_000_000L) {
				flush();
			}
		}

		void flush() {

			try (var out = new FileWriter(file, true)) {

				for (var entry : buf) {
					out.write(entry);
				}
				buf.clear();
				failed = false;

			} catch (IOException ex) {
				if (!failed) {
					ex.printStackTrace(System.err);
					failed = true;
				}
			}

			flushNs = System.nanoTime();
		}
	}

	public final MultiStateConfSpace confSpace;

	private final StatePerf[] states;

	private PerfLog perfLog = null;

	public NodePerformance(MultiStateConfSpace confSpace) {

		this.confSpace = confSpace;

		states = confSpace.states.stream()
			.map(state -> new StatePerf(state))
			.toArray(StatePerf[]::new);
	}

	public synchronized void clear() {
		for (var state : states) {
			state.clear();
		}
	}

	public void setLog(File file) {
		if (file != null) {
			perfLog = new PerfLog(file);
		} else {
			perfLog = null;
		}
	}

	public synchronized void updateAndLog(int statei, int[] conf, long ns, BigExp zSumUpper, BigExp reduction, BigExp score) {
	    updateAndLog(statei, conf, ns, zSumUpper, reduction, score, new BigExp(0.0), new BigExp(0.0), new BigExp(0.0));
    }

	public synchronized void updateAndLog(int statei, int[] conf, long ns, BigExp zSumUpper, BigExp reduction, BigExp score, BigExp zSumLower, BigExp zSumUpperFactor, BigExp zSumLowerFactor) {

		states[statei].get(conf).update(ns);

		if (perfLog != null) {
			perfLog.add(statei, conf, ns, zSumUpper, reduction, score, zSumLower, zSumUpperFactor, zSumLowerFactor);
		}
	}

	public synchronized void update(int statei, int[] conf, long ns) {
		states[statei].get(conf).update(ns);
	}

	public void updateAndLog(NodeIndex.Node node, long ns, BigExp reduction) {
		updateAndLog(node.statei, node.conf, ns, node.zSumUpper, reduction, node.score, node.zSumLower, node.zSumUpperFactor, node.zSumLowerFactor);
	}

	public synchronized BigExp score(int statei, int[] conf, BigExp zSumUpper) {
		return states[statei].get(conf).score(zSumUpper);
	}

	public BigExp score(NodeIndex.Node node) {
		return score(node.statei, node.conf, node.zSumUpper);
	}
}
