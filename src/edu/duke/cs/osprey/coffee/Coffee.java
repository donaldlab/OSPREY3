package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.BigExp;

import java.math.MathContext;
import java.util.Arrays;
import java.util.function.Consumer;


/**
 * A clusterized successor to SOFEA that uses compiled conformation spaces.
 *
 * It's basically just a fancy implementation of MARK*
 * that doesn't rely on the traditional DEE->A* pipeline
 */
public class Coffee {

	public static class Builder {

		public final MultiStateConfSpace confSpace;

		private final StateConfig[] stateConfigs;
		private Cluster cluster;
		private Parallelism parallelism;

		public Builder(MultiStateConfSpace confSpace) {
			this.confSpace = confSpace;
			stateConfigs = new StateConfig[confSpace.states.size()];
		}

		public Builder configState(MultiStateConfSpace.State state, StateConfig config) {
			stateConfigs[state.index] = config;
			return this;
		}

		public Builder configEachState(Consumer<StateConfig> configurator) {
			for (MultiStateConfSpace.State state : confSpace.states) {
				StateConfig config = new StateConfig(state);
				configurator.accept(config);
				configState(state, config);
			}
			return this;
		}

		public Builder setCluster(Cluster val) {
			cluster = val;
			return this;
		}

		public Builder setParallelism(Parallelism val) {
			parallelism = val;
			return this;
		}

		public Coffee build() {

			// check the state configs
			for (var config : stateConfigs) {
				config.check();
			}

			// make a single-node cluster if needed
			if (cluster == null) {
				cluster = new Cluster("COFFEE", "job", 0, 1);
			}

			// make default parallelism if needed
			if (parallelism == null) {
				parallelism = Parallelism.makeCpu(1);
			}

			return new Coffee(confSpace, stateConfigs, cluster, parallelism);
		}
	}


	public static class StateConfig {

		public final MultiStateConfSpace.State state;

		public ConfEnergyCalculator ecalc;
		public PosInterGen posInterGen;

		public StateConfig(MultiStateConfSpace.State state) {
			this.state = state;
		}

		public void check() {
			if (ecalc == null) {
				throw new MissingArgumentException("conformation energy calculator (ecalc)", state);
			}
			if (posInterGen == null) {
				throw new MissingArgumentException("position interactions generator (posInterGen)", state);
			}
		}

		static class MissingArgumentException extends RuntimeException {
			MissingArgumentException(String name, MultiStateConfSpace.State state) {
				super("Please set a value for: " + name + " for the state: " + state.name);
			}
		}
	}


	/**
	 * Tells when to stop refining
	 */
	public interface Criterion {
		boolean isFinished();
	}

	public final MultiStateConfSpace confSpace;
	public final StateConfig[] stateConfigs;
	public final Cluster cluster;
	public final Parallelism parallelism;

	public final StateInfo[] infos;

	public final MathContext mathContext = BigExp.mathContext;
	public final BoltzmannCalculator bcalc = new BoltzmannCalculator(mathContext);

	private Coffee(MultiStateConfSpace confSpace, StateConfig[] stateConfigs, Cluster cluster, Parallelism parallelism) {
		this.confSpace = confSpace;
		this.stateConfigs = stateConfigs;
		this.cluster = cluster;
		this.parallelism = parallelism;
		infos = Arrays.stream(stateConfigs)
			.map(config -> new StateInfo(config, bcalc))
			.toArray(StateInfo[]::new);
	}

	public void refine() {
		try (var member = new ClusterMember(cluster)) {
			try (var tasks = parallelism.makeTaskExecutor()) {

				NodeDB nodedb = new NodeDB(confSpace, member, null);

				// wait for everyone to get here
				member.log0("waiting for cluster to assemble ...");
				member.barrier();

				// pre-compute the Z matrices
				for (var info : infos) {
					member.log0("computing Z matrix for state: %s", info.config.state.name);
					info.zmat.compute(member, tasks);
				}

				// compute bounds on free energies at the root nodes
				for (int statei=0; statei<confSpace.states.size(); statei++) {
					var stateInfo = infos[statei];

					// get a multi-sequence Z bound on the root node
					ConfIndex index = stateInfo.makeConfIndex();
					BigExp zSumUpper = stateInfo.zSumUpper(index);

					member.log0("state: %10s, zSumUpper: %s", confSpace.states.get(statei).name, zSumUpper);

					// make sure BigExp values are fully normalized before writing to the databases to avoid some roundoff error
					zSumUpper.normalize(true);

					// init the nodedb with the root node
					nodedb.init(statei, Conf.make(index), zSumUpper);

					// TODO: init sequence database?
					//seqtx.addZSumUpper(state, state.confSpace.seqSpace().makeUnassignedSequence(), zSumUpper);
				}

				/* TODO: criterion, sweep thresholds, and refinement
				// init sweep thresholds for each state
				Double[] gThresholds = confSpace.states.stream()
					.map(state -> {
						return 5.0; // TODO
					})
					.toArray(Double[]::new);

				while (true) {

					// check the criterion
					if (criterion.isFinished()) {
						break;
					}

					// update the sweep threshold
					for (MultiStateConfSpace.State state : confSpace.states) {
						if (gThresholds[state.index] != null) {
							if (gThresholds[state.index] > gThresholdUpper) {
								gThresholds[state.index] = null;
							} else {
								gThresholds[state.index] += sweepIncrement;
							}
						}
					}
				}
				*/
			}
		}
	}
}
