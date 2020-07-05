package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.coffee.nodedb.NodeDB;
import edu.duke.cs.osprey.coffee.nodedb.NodeIndex;
import edu.duke.cs.osprey.coffee.seqdb.SeqDB;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.gpu.Structs;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.parallelism.ThreadTools;
import edu.duke.cs.osprey.tools.BigExp;

import java.io.File;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.TimeUnit;
import java.util.function.BiConsumer;


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
		private Structs.Precision precision = Structs.Precision.Float64;
		private File nodedbFile = null;
		private long nodedbFileBytes = 0;
		private long nodedbMemBytes = 2*1024*1024; // 2 MiB
		private File seqdbFile = null;
		private MathContext seqdbMathContext = new MathContext(128, RoundingMode.HALF_UP);
		private boolean includeStaticStatic = true;
		private BoltzmannCalculator.Conditions conditions = BoltzmannCalculator.Conditions.Classic; // don't rock the boat
		// TODO: experiment if other conditions better correlate with experimental results?

		public Builder(MultiStateConfSpace confSpace) {
			this.confSpace = confSpace;
			stateConfigs = new StateConfig[confSpace.states.size()];
		}

		public Builder configState(MultiStateConfSpace.State state, BiConsumer<StateConfig,ConfEnergyCalculator> configurator) {
			StateConfig config = new StateConfig(state);
			var ecalc = new CPUConfEnergyCalculator(config.confSpace);
			configurator.accept(config, ecalc);
			stateConfigs[state.index] = config;
			return this;
		}

		public Builder configState(String stateName, BiConsumer<StateConfig,ConfEnergyCalculator> configurator) {
			return configState(confSpace.getState(stateName), configurator);
		}

		public Builder configEachState(BiConsumer<StateConfig,ConfEnergyCalculator> configurator) {
			for (MultiStateConfSpace.State state : confSpace.states) {
				configState(state, configurator);
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

		public Builder setPrecision(Structs.Precision val) {
			precision = val;
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

		public Builder setStaticStatic(boolean val) {
			includeStaticStatic = val;
			return this;
		}

		public Builder setConditions(BoltzmannCalculator.Conditions val) {
			conditions = val;
			return this;
		}

		public Builder setTemp(double val) {
			conditions = new BoltzmannCalculator.Conditions(val);
			return this;
		}

		public Coffee build() {

			// check the state configs
			for (var config : stateConfigs) {
				config.check();
			}

			// make a single-node cluster if needed
			if (cluster == null) {
				// if no cluster was given, make a random id,
				// so we (probably) don't try to join with other jobs on this machine, and vice versa
				int jobId = Math.abs(new Random().nextInt());
				cluster = new Cluster("COFFEE", "job-" + jobId, 0, 1);
			}

			// make default parallelism if needed
			if (parallelism == null) {
				parallelism = Parallelism.makeCpu(1);
			}

			return new Coffee(
				confSpace, stateConfigs, cluster, parallelism, precision,
				nodedbFile, nodedbFileBytes, nodedbMemBytes,
				seqdbFile, seqdbMathContext, includeStaticStatic, conditions
			);
		}
	}


	public static class StateConfig {

		public final MultiStateConfSpace.State state;
		public final ConfSpace confSpace;

		public PosInterGen posInterGen;

		public StateConfig(MultiStateConfSpace.State state) {

			if (!(state.confSpace instanceof ConfSpace)) {
				throw new IllegalArgumentException("Conformation space must be in the compiled format");
			}

			this.state = state;
			this.confSpace = (ConfSpace)state.confSpace;
		}

		public void check() {
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

		/**
		 * Return the number of best conformations to track for each state and sequence.
		 */
		default int numBestConfs() {
			return 0;
		}

		/**
		 * Prepare all the databases and cluster members for computation.
		 */
		void init(Directions directions, NodeProcessor processor);

		/**
		 * Start the computation, manage the cluster members, and finish the computation.
		 */
		void direct(Directions directions, NodeProcessor processor);
	}

	public final MultiStateConfSpace confSpace;
	public final StateConfig[] stateConfigs;
	public final Cluster cluster;
	public final Parallelism parallelism;
	public final Structs.Precision precision;
	public final File dbFile;
	public final long dbFileBytes;
	public final long dbMemBytes;
	public final File seqdbFile;
	public final MathContext seqdbMathContext;
	public final boolean includeStaticStatic;
	public final BoltzmannCalculator.Conditions conditions;

	public final MathContext mathContext = BigExp.mathContext;
	public final BoltzmannCalculator bcalc;
	public final StateInfo[] infos;

	private Coffee(MultiStateConfSpace confSpace, StateConfig[] stateConfigs, Cluster cluster, Parallelism parallelism, Structs.Precision precision, File dbFile, long dbFileBytes, long dbMemBytes, File seqdbFile, MathContext seqdbMathContext, boolean includeStaticStatic, BoltzmannCalculator.Conditions conditions) {

		this.confSpace = confSpace;
		this.stateConfigs = stateConfigs;
		this.cluster = cluster;
		this.parallelism = parallelism;
		this.precision = precision;
		this.dbFile = dbFile;
		this.dbFileBytes = dbFileBytes;
		this.dbMemBytes = dbMemBytes;
		this.seqdbFile = seqdbFile;
		this.seqdbMathContext = seqdbMathContext;
		this.includeStaticStatic = includeStaticStatic;
		this.conditions = conditions;

		bcalc = new BoltzmannCalculator(mathContext, conditions);
		infos = Arrays.stream(stateConfigs)
			.map(config -> new StateInfo(config, bcalc))
			.toArray(StateInfo[]::new);
	}

	public void run(Director director) {
		try (var member = new ClusterMember(cluster)) {

			// wait for everyone to get here
			member.log0("waiting for cluster to assemble ...");
			member.barrier(1, TimeUnit.MINUTES);

			try (var tasks = new ThreadPoolTaskExecutor()) {
				tasks.start(parallelism.numThreads);

				// open the sequence database
				try (SeqDB seqdb = new SeqDB.Builder(confSpace, member)
					.setMathContext(seqdbMathContext)
					.setFile(seqdbFile)
					.setNumBestConfs(director.numBestConfs())
					.build()
				) {

					// open the node database
					try (var nodedb = new NodeDB.Builder(confSpace, member)
						.setFile(dbFile, dbFileBytes)
						.setMem(dbMemBytes)
						.build()
					) {

						// init the node processor, and report dropped nodes to the sequence database
						try (var nodeProcessor = new NodeProcessor(tasks, seqdb, nodedb, infos, includeStaticStatic, parallelism, precision)) {
							nodedb.dropHandler = nodeProcessor::handleDrops;

							// wait for everyone to be ready
							member.barrier(5, TimeUnit.MINUTES);

							// pre-compute the Z matrices
							for (var info : infos) {
								member.log0("computing Z matrix for state: %s", info.config.state.name);
								info.zmat.compute(member, tasks, includeStaticStatic, nodeProcessor.ecalcs[info.config.state.index]);
							}

							// initialize the directions and wait
							var directions = new Directions(member);
							member.barrier(5, TimeUnit.MINUTES);

							// initialize the computation
							if (member.isDirector()) {
								director.init(directions, nodeProcessor);
								initRootNodes(directions, nodeProcessor);
							}

							// wait for the root nodes calculation to finish
							// TODO: do we need to configure the init time here for different init procedures?
							member.barrier(5, TimeUnit.MINUTES);

							// prep complete! now we can start the real computation
							if (member.isDirector()) {
								director.direct(directions, nodeProcessor);
							} else {
								followDirections(directions, nodeProcessor);
							}

							// wait for the computation to finish before cleaning up databases
							member.barrier(5, TimeUnit.MINUTES);

						} // node processor
					} // nodedb
				} // seqdb
			} // tasks
		} // cluster member
	}

	private void initRootNodes(Directions directions, NodeProcessor processor) {

		// compute bounds on free energies at the root nodes
		var batch = processor.seqdb.batch();
		for (int statei=0; statei<confSpace.states.size(); statei++) {
			var stateInfo = infos[statei];

			// get a (possibly) multi-sequence Z bound on the root node
			ConfIndex index = stateInfo.makeConfIndex();
			var tree = directions.getTreeOrThrow(statei);
			BigExp zSumUpper = stateInfo.zSumUpper(index, tree).normalize(true);

			// init the node database
			var node = new NodeIndex.Node(statei, Conf.make(index), zSumUpper, zSumUpper);
			processor.nodedb.addLocal(node);

			// init sequence database
			batch.addZSumUpper(
				stateInfo.config.state,
				stateInfo.config.state.isSequenced ? confSpace.seqSpace.makeUnassignedSequence() : null,
				node.score
			);
		}
		batch.save();

		// let the rest of the cluster know right away we have the root nodes
		processor.nodedb.broadcast();
	}

	private void followDirections(Directions directions, NodeProcessor processor) {

		while (directions.isRunning()) {

			// try to process a node
			// TEMP
			//var result = processor.process(directions);
			var result = NodeProcessor.Result.NoState;

			switch (result) {

				// TODO

				// wait a bit and try again
				case NoState -> ThreadTools.sleep(500, TimeUnit.MILLISECONDS);
			}
		}
	}
}
