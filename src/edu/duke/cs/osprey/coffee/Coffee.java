package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.coffee.commands.Commands;
import edu.duke.cs.osprey.coffee.nodedb.NodeDB;
import edu.duke.cs.osprey.coffee.nodedb.NodeIndex;
import edu.duke.cs.osprey.coffee.seqdb.SeqDB;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.BigExp;

import java.io.File;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;


/**
 * A clusterized successor to SOFEA that uses compiled conformation spaces.
 *
 * It's basically just a fancy implementation of MARK*
 * that has better memory usage than the traditional DEE->A* pipeline
 */
public class Coffee {

	public static class Builder {

		public final MultiStateConfSpace confSpace;

		private final StateConfig[] stateConfigs;
		private Cluster cluster;
		private Parallelism parallelism;
		private File nodedbFile = null;
		private long nodedbFileBytes = 0;
		private long nodedbMemBytes = 2*1024*1024; // 2 MiB
		private File seqdbFile = null;
		private MathContext seqdbMathContext = new MathContext(128, RoundingMode.HALF_UP);

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

		public Builder setNodeDBFile(File file, long bytes) {
			nodedbFile = file;
			nodedbFileBytes = bytes;
			return this;
		}

		public Builder setNodeDBMem(long bytes) {
			nodedbMemBytes = bytes;
			return this;
		}

		public Builder setSeqDBFile(File file) {
			seqdbFile = file;
			return this;
		}

		public Builder setSeqDBMathContext(MathContext val) {
			seqdbMathContext = val;
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

			return new Coffee(confSpace, stateConfigs, cluster, parallelism, nodedbFile, nodedbFileBytes, nodedbMemBytes, seqdbFile, seqdbMathContext);
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
	 * Tells COFFEE what nodes to look for and when to stop looking
	 */
	public interface Director {
		void direct(Commands commands, NodeProcessor processor);
	}

	public final MultiStateConfSpace confSpace;
	public final StateConfig[] stateConfigs;
	public final Cluster cluster;
	public final Parallelism parallelism;
	public final File dbFile;
	public final long dbFileBytes;
	public final long dbMemBytes;
	public final File seqdbFile;
	public final MathContext seqdbMathContext;

	public final StateInfo[] infos;

	public final MathContext mathContext = BigExp.mathContext;
	public final BoltzmannCalculator bcalc = new BoltzmannCalculator(mathContext);

	private Coffee(MultiStateConfSpace confSpace, StateConfig[] stateConfigs, Cluster cluster, Parallelism parallelism, File dbFile, long dbFileBytes, long dbMemBytes, File seqdbFile, MathContext seqdbMathContext) {
		this.confSpace = confSpace;
		this.stateConfigs = stateConfigs;
		this.cluster = cluster;
		this.parallelism = parallelism;
		this.dbFile = dbFile;
		this.dbFileBytes = dbFileBytes;
		this.dbMemBytes = dbMemBytes;
		this.seqdbFile = seqdbFile;
		this.seqdbMathContext = seqdbMathContext;
		infos = Arrays.stream(stateConfigs)
			.map(config -> new StateInfo(config, bcalc))
			.toArray(StateInfo[]::new);
	}

	public void run(Director director) {
		try (var member = new ClusterMember(cluster)) {
			try (var tasks = parallelism.makeTaskExecutor()) {

				// wait for everyone to get here
				member.log0("waiting for cluster to assemble ...");
				member.barrier(1, TimeUnit.MINUTES);

				// pre-compute the Z matrices
				for (var info : infos) {
					member.log0("computing Z matrix for state: %s", info.config.state.name);
					info.zmat.compute(member, tasks);
				}

				// open the sequence database
				try (SeqDB seqdb = new SeqDB(confSpace, seqdbMathContext, seqdbFile, member)) {

					// open the node database
					try (var nodedb = new NodeDB.Builder(confSpace, member)
						.setFile(dbFile, dbFileBytes)
						.setMem(dbMemBytes)
						.build()
					) {

						if (member.isDirector()) {

							var batch = seqdb.batch();

							// compute bounds on free energies at the root nodes
							for (int statei=0; statei<confSpace.states.size(); statei++) {
								var stateInfo = infos[statei];

								// start with the static-static energy
								BigExp zSumUpper = stateInfo.zmat.staticStatic();

								// get a multi-sequence Z bound on the root node
								ConfIndex index = stateInfo.makeConfIndex();
								zSumUpper.mult(stateInfo.zSumUpper(index));

								// make sure BigExp values are fully normalized before writing to the databases to avoid some roundoff error
								zSumUpper.normalize(true);

								// init the nodedb with the root node
								nodedb.addLocal(new NodeIndex.Node(statei, Conf.make(index), zSumUpper));

								// init sequence database
								batch.addZSumUpper(
									stateInfo.config.state,
									confSpace.seqSpace.makeUnassignedSequence(),
									zSumUpper
								);
							}

							batch.save();

							// TEMP
							member.log("seqdb: %s", seqdb.dump());

							// let the rest of the cluster know we have the root nodes
							nodedb.broadcast();
						}

						// wait for the root nodes calculation to finish
						member.barrier(1, TimeUnit.MINUTES);

						// TEMP
						member.log("starting computation");

						// prep complete! now we can start the real computation
						var commands = new Commands(member);
						var nodeProcessor = new NodeProcessor(tasks, seqdb, nodedb);
						if (member.isDirector()) {
							director.direct(commands, nodeProcessor);
						} else {
							follow(commands, nodeProcessor);
						}

						// TEMP
						member.log("finished computation");

						// wait for the computation to finish before cleaning up databases
						member.barrier(5, TimeUnit.MINUTES);

					} // nodedb
				} // seqdb
			} // tasks
		} // cluster member
	}

	private void follow(Commands commands, NodeProcessor processor) {

		while (commands.isRunning()) {

			// get the currently focused state, or wait for one to happen
			int statei = commands.getFocusedStatei();
			if (statei < 0) {
				commands.member.sleep(500, TimeUnit.MILLISECONDS);
				continue;
			}

			// get the next node from that staet, or wait for one to happen
			NodeIndex.Node node = processor.nodedb.removeHigh(statei);
			if (node == null) {
				commands.member.sleep(500, TimeUnit.MILLISECONDS);
				continue;
			}

			// TODO: use tasks?
			processor.process(node);
		}
	}
}
