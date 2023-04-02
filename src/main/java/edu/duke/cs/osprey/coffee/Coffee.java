package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.coffee.nodedb.NodeDB;
import edu.duke.cs.osprey.coffee.nodedb.NodeIndex;
import edu.duke.cs.osprey.coffee.nodedb.NodeTree;
import edu.duke.cs.osprey.coffee.seqdb.SeqDB;
import edu.duke.cs.osprey.coffee.zmat.ClusterZMatrix;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.NativeConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.gpu.Structs;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.math.MathContext;
import java.math.RoundingMode;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
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

		/** Information about the cluter resources available, if any. */
		private Cluster cluster;

		/** Information about parallel hardware resources available, if any. */
		private Parallelism parallelism;

		/** The precision of floating-point calculations for the energy function. */
		private Structs.Precision precision = Structs.Precision.Float64;

		/** Path to the Node Database file, if any. If no file is given, the calculation will be done in RAM. */
		private File nodedbFile = null;

		/** Size of the Node Database file. */
		private long nodedbFileBytes = 0;

		/** Amount of memory (RAM, in bytes) to use for the calculation. This memory will be pre-allocated. */
		private long nodedbMemBytes = 2*1024*1024; // 2 MiB

		/** Path to the file for the Sequence Database, if any. If no file is given, the sequence info will be stored in RAM. */
		private File seqdbFile = null;

		/**
		 * Precision to use for calculations involving sequence partition function values.
		 *
		 * @note SeqDB needs a LOT of precision to keep the bounds acccurate,
		 * since we're adding/subtracting numbers with very different exponents
		 * empirical testing shows anything higher than 512 starts to show noticeable performance penalties
		 * unfortunately, testing also shows that values of 1024 of less can cause some pfunc computations to get stuck
		 * so we'll need something really huge by default to be safe for most computations
		 */
		private MathContext seqdbMathContext = new MathContext(2048, RoundingMode.HALF_UP);

		/**
		 * True to include interactions within the static region (unaffected by motions and mutations) of the input molecules.
		 *
		 * Including these energies, or not, won't change the rankings of sequences/conformations in a single design,
		 * but including them will make your energies comaparable across different designs on the same molecules.
		 *
		 * Osprey has been optimized so that including these interactions, or not, should have no noticeable
		 * impact on performance.
		 */
		private boolean includeStaticStatic = true;

		/**
		 * Threshold (in kcal/mol) for when triple interactions (ie between three design positions, rather than one or two) should
		 * be included in the energy matrix.
		 *
		 * Use null to never include triple corrections to the energy matrix.
		 *
		 * Use a non-null value to include triple corrections in the energy matrix.
		 * For example, a value like 10 kcal/mol seems to work well, but feel free to experiment with different values
		 * in your own designs.
		 *
		 * Pairwise energies in energy matrices are often overly optimistic (ie too low) which leads to loose
		 * lower bounds on conformation energies. Using triple information to correct (ie raise) the energy
		 * bounds on conformations can speed up large Osprey calculations significantly by using a more accurate
		 * conformtion order to compute free energies. However, the tradeoff for faster free energy calculation is
		 * a slower energy matrix calculation step, so the speedup tends to only pay off for larger designs.
		 *
		 * A triple correction will be added to the energy matix when all of the triple's constituent single and pairwise energies
		 * are below the threshold. This setting ensures we don't waste time correcting energy matrix energies that are
		 * already very high, like clashes, since triple corrections are expensive to compute.
		 */
		private Double tripleCorrectionThreshold = null;

		/**
		 * Environmental conditions (like temperature) for calculations that use them (like Boltzmann weighting).
		 */
		private BoltzmannCalculator.Conditions conditions = BoltzmannCalculator.Conditions.Classic; // don't rock the boat
		// TODO: experiment if other conditions better correlate with experimental results?

		/**
		 * For debugging Osprey performance, most users will not need this.
		 */
		private File nodeScoringLog = null;

		/**
		 * For debugging Osprey performance, most users will not need this.
		 */
		private Duration nodeStatsReportingInterval = null;

		public Builder(MultiStateConfSpace confSpace) {
			this.confSpace = confSpace;
			stateConfigs = new StateConfig[confSpace.states.size()];
		}

		public Builder configState(StateConfig config) {
			stateConfigs[config.state.index] = config;
			return this;
		}

		public Builder configState(MultiStateConfSpace.State state, Consumer<StateConfig> configurator) {
			StateConfig config = new StateConfig(state);
			configurator.accept(config);
			return configState(config);
		}

		public Builder configState(String stateName, Consumer<StateConfig> configurator) {
			return configState(confSpace.getState(stateName), configurator);
		}

		public Builder configEachState(Consumer<StateConfig> configurator) {
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

		public Builder setTripleCorrectionThreshold(Double val) {
			tripleCorrectionThreshold = val;
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

		public Builder setNodeScoringLog(File val) {
			nodeScoringLog = val;
			return this;
		}

		public Builder setNodeStatsReportingInterval(Duration val) {
			nodeStatsReportingInterval = val;
			return this;
		}

		public Coffee build() {

			// check the state configs
			for (var config : stateConfigs) {
				config.check();
			}

			// try to auto-detect the cluster
			if (cluster == null) {
				cluster = Cluster.fromSLURM();
			}
			if (cluster == null) {

				// if no cluster was given/found, make a random id,
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
				seqdbFile, seqdbMathContext, includeStaticStatic, tripleCorrectionThreshold,
				conditions, nodeScoringLog, nodeStatsReportingInterval
			);
		}
	}


	public static class StateConfig {

		public final MultiStateConfSpace.State state;
		public final ConfSpace confSpace;

		// provide a simple default that should be fine for most cases
		public PosInterGen posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, null);

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

		public List<PosInter> makeInters(int[] conf, boolean includeStaticStatic) {
			if (includeStaticStatic) {
				return posInterGen.all(confSpace, conf);
			} else {
				return posInterGen.dynamic(confSpace, conf);
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
	public final Double tripleCorrectionThreshold;
	public final BoltzmannCalculator.Conditions conditions;
	public final File nodeScoringLog;
	public final Duration nodeStatsReportingInterval;

	public final MathContext mathContext = BigExp.mathContext;
	public final BoltzmannCalculator bcalc;
	public final StateInfo[] infos;

	private Coffee(
		MultiStateConfSpace confSpace, StateConfig[] stateConfigs, Cluster cluster, Parallelism parallelism, Structs.Precision precision,
		File dbFile, long dbFileBytes, long dbMemBytes,
		File seqdbFile, MathContext seqdbMathContext, boolean includeStaticStatic, Double tripleCorrectionThreshold,
		BoltzmannCalculator.Conditions conditions, File nodeScoringLog, Duration nodeStatsReportingInterval
	) {

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
		this.tripleCorrectionThreshold = tripleCorrectionThreshold;
		this.conditions = conditions;
		this.nodeScoringLog = nodeScoringLog;
		this.nodeStatsReportingInterval = nodeStatsReportingInterval;

		bcalc = new BoltzmannCalculator(mathContext, conditions);
		infos = Arrays.stream(stateConfigs)
			.map(config -> new StateInfo(config))
			.toArray(StateInfo[]::new);
	}

	public void run(Director director) {
		try (var member = new ClusterMember(cluster)) {

			var stopwatch = new Stopwatch().start();

			// wait for everyone to get here
			member.log0("waiting for cluster to assemble ...");
			member.barrier(1, TimeUnit.MINUTES);

			try (var cpuTasks = new ThreadPoolTaskExecutor()) {
				cpuTasks.start(parallelism.numThreads);

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
						.setScoringLog(nodeScoringLog)
						.build()
					) {

						// init the node processor, and report dropped nodes to the sequence database
						try (var nodeProcessor = new NodeProcessor(cpuTasks, seqdb, nodedb, infos, includeStaticStatic, parallelism, precision, nodeStatsReportingInterval)) {
							nodedb.setDropHandler(nodeProcessor::handleDrops);

							// wait for everyone to be ready
							member.barrier(5, TimeUnit.MINUTES);

							// pre-compute the Z matrices
							for (var info : infos) {
								member.log0("computing Z matrix for state: %s", info.config.state.name);
								ConfEnergyCalculator ecalc = nodeProcessor.cpuEcalcs[info.config.state.index];
								if (nodeProcessor.gpuEcalcs != null) {
									ecalc = nodeProcessor.gpuEcalcs[info.config.state.index];
								}
								info.zmat = new ClusterZMatrix(info.config.confSpace, info.config.posInterGen, bcalc);
								info.zmat.compute(member, cpuTasks, includeStaticStatic, tripleCorrectionThreshold, ecalc);
								info.initBounder();
							}

							// initialize the directions and wait
							var directions = new Directions(confSpace, member);
							member.barrier(5, TimeUnit.MINUTES);

							// start the node threads
							nodeProcessor.start(parallelism.numThreads, directions);
							try {

								// prep complete! now we can start the real computation
								if (member.isDirector()) {
									director.direct(directions, nodeProcessor);
								} else {
									directions.waitForStop();
								}

								member.log0("Computation finished, cleaning up resources ...");

							} catch (Throwable t) {

								// catch and report errors before calling close() on the node processor
								// since calling close() without cleaning up gracefully can cause deadlock

								t.printStackTrace(System.err);

								member.log("Error encountered, stopping cluster ...");
							}

							// ask the cluster to stop nicely
							directions.stop();

							// wait for the computation to finish before cleaning up databases
							member.barrier(5, TimeUnit.MINUTES);

						} // node processor
					} // nodedb
				} // seqdb
			} // tasks

			member.log0("COFFEE Finished in %s", stopwatch.getTime(2));

		} // cluster member
	}

	/**
	 * Compute a Z matrix.
	 * Mostly only useful for debugging.
	 */
	public ClusterZMatrix calcZMat(int statei) {

		var stateInfo = infos[statei];

		try (var member = new ClusterMember(cluster)) {
			try (var cpuTasks = new ThreadPoolTaskExecutor()) {
				cpuTasks.start(parallelism.numThreads);

				try (var ecalc = new CPUConfEnergyCalculator(stateInfo.config.confSpace)) {

					var zmat = new ClusterZMatrix(stateInfo.config.confSpace, stateInfo.config.posInterGen, bcalc);
					zmat.compute(member, cpuTasks, includeStaticStatic, tripleCorrectionThreshold, ecalc);
					return zmat;
				}
			}
		}
	}

	/**
	 * Quickly get a few nodes with high Z values.
	 * Mostly only useful for debugging.
	 */
	public List<NodeIndex.Node> findHighZNodes(int statei, ClusterZMatrix zmat, NodeTree tree, int count) {

		var leafNodes = new ArrayList<NodeIndex.Node>(count);

		int batchSize = 10;
		var nodeBatch = new ArrayList<NodeIndex.Node>(batchSize);
		var expandedNodes = new ArrayList<NodeIndex.Node>();

		try (var member = new ClusterMember(cluster)) {

			try (var cpuTasks = new ThreadPoolTaskExecutor()) {
				cpuTasks.start(parallelism.numThreads);

				// open the node database
				try (var nodedb = new NodeDB.Builder(confSpace, member)
					.setFile(dbFile, dbFileBytes)
					.setMem(dbMemBytes)
					.setScoringLog(nodeScoringLog)
					.build()
				) {

					// init the node processor, and report dropped nodes to the sequence database
					try (var nodeProcessor = new NodeProcessor(cpuTasks, null, nodedb, infos, includeStaticStatic, parallelism, precision, nodeStatsReportingInterval)) {

						// init the state with the zmat
						var stateInfo = infos[statei];
						stateInfo.zmat = zmat;
						stateInfo.initBounder();

						// init the root node
						ConfIndex index = stateInfo.makeConfIndex();
						BigExp zSumUpper = stateInfo.zSumUpper(index, tree).normalize(true);
						nodeProcessor.nodedb.addLocal(new NodeIndex.Node(statei, Conf.make(index), zSumUpper, zSumUpper));

						while (true) {

							// get the next nodes
							nodedb.removeHighestLocal(statei, batchSize, nodeBatch);
							if (nodeBatch.isEmpty()) {
								throw new RuntimeException("couldn't find " + count + " nodes");
							}
							for (var node : nodeBatch) {

								// if it's a leaf node, add it to the list
								if (node.isLeaf()) {
									leafNodes.add(node);
									if (leafNodes.size() >= count) {
										return leafNodes;
									}
									continue;
								}

								// otherwise, expand the node
								nodeProcessor.expand(node, tree, expandedNodes);
							}
							nodeBatch.clear();

							// add the nodes to the db
							nodedb.addLocal(statei, expandedNodes);
							expandedNodes.clear();
						}

					} // node processor
				} // nodedb
			} // cpu tasks
		} // cluster member
	}
}
