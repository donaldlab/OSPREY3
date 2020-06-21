package edu.duke.cs.osprey.coffee.nodedb;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.BigExp;

import java.util.ArrayDeque;
import java.util.Deque;


public class NodePerformance {

	public static final long InitialNs = 100_000L; // start with 100 us per node operation
	public static final int HistorySize = 100;

	private static final BigExp zero = new BigExp(0.0, 0);
	private static final BigExp one = new BigExp(1.0, 0);

	private static class TypePerf {

		Deque<Long> nsHistory = new ArrayDeque<>(HistorySize);
		Deque<BigExp> reductionFactorHistory = new ArrayDeque<>(HistorySize);
		long nsSum = 0L;
		BigExp reductionFactorSum = new BigExp(zero);
		long avgNs = InitialNs;
		BigExp avgReductionFactor = new BigExp(one);

		TypePerf() {

			// init the history with something neutral, so we can ease-in changes over time
			for (int i=0; i<HistorySize; i++) {
				nsHistory.add(InitialNs);
				nsSum += InitialNs;
				reductionFactorHistory.add(one);
				reductionFactorSum.add(one);
			}
		}

		void update(long ns, BigExp zSumUpper, BigExp reduction) {

			// add the new time
			nsSum -= nsHistory.removeFirst();
			nsHistory.add(ns);
			nsSum += ns;

			// add the new reduction factor
			var reductionFactor = new BigExp(reduction);
			reductionFactor.div(zSumUpper);
			reductionFactorSum.sub(reductionFactorHistory.removeFirst());
			reductionFactorHistory.add(reductionFactor);
			reductionFactorSum.add(reductionFactor);

			// update the averages
			avgNs = nsSum/HistorySize;
			avgReductionFactor.set(reductionFactorSum);
			avgReductionFactor.div(HistorySize);
		}

		public BigExp predictReduction(BigExp zSumUpper) {
			var reduction = new BigExp(zSumUpper);
			reduction.mult(avgReductionFactor);
			return reduction;
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

		TypePerf get(int[] conf) {
			return types[Conf.countAssignments(conf)];
		}
	}

	public final MultiStateConfSpace confSpace;

	private final StatePerf[] states;

	public NodePerformance(MultiStateConfSpace confSpace) {

		this.confSpace = confSpace;

		states = confSpace.states.stream()
			.map(state -> new StatePerf(state))
			.toArray(StatePerf[]::new);
	}

	public synchronized void update(int statei, int[] conf, long ns, BigExp zSumUpper, BigExp reduction) {
		states[statei].get(conf).update(ns, zSumUpper, reduction);
	}

	public void update(NodeIndex.Node node, long ns, BigExp reduction) {
		update(node.statei, node.conf, ns, node.zSumUpper, reduction);
	}

	public synchronized BigExp score(int statei, int[] conf, BigExp zSumUpper) {
		var type = states[statei].get(conf);
		var score = type.predictReduction(zSumUpper);
		score.div(type.avgNs);
		score.normalize(true);
		return score;
	}

	public BigExp score(NodeIndex.Node node) {
		return score(node.statei, node.conf, node.zSumUpper);
	}
}
