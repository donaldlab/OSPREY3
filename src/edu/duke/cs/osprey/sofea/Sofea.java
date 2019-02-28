package edu.duke.cs.osprey.sofea;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.*;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;


/**
 * SOFEA - Sweep Operations for Free Energy Approximation
 *
 * Explores a (possibly multi-sequence and multi-state) conformation space by visiting
 * subtrees containing low-energy conformations first, in less than exponential memory.
 * Keeps a running free-energy calculation going for each state and sequence.
 * Free energy approximations are improved by sweeping an energy threshold over a bounded
 * parameter space. If the entire parameter space is swept, the free energies are guaranteed
 * to be completely precise, but hopefully you can stop before then.
 *
 * SOFEA can be used to do multi-state design when used with a stopping criterion that stops
 * the sweep when enough sequences are found whose free energies optimize some objective function.
 */
public class Sofea {

	public static class Builder {

		public final MultiStateConfSpace confSpace;

		private StateConfig[] stateConfigs;

		/**
		 * SOFEA does math with very large floating numbers, how much decimal precision should we use?
		 */
		private MathContext mathContext = new MathContext(32, RoundingMode.HALF_UP);

		/**
		 * File for the sequence database
		 */
		private File seqdbFile = new File("seq.db");

		/**
		 * The sequence database adds very large and very small numbers together,
		 * so it needs a bit more floating point precision.
		 */
		private MathContext seqdbMathContext = new MathContext(128, RoundingMode.HALF_UP);

		/**
		 * File for the fringe datbase, ie the set of fringe nodes
		 */
		private File fringedbFile = new File("fringe.db");

		/**
		 * Max size of the fringe database, in bytes
		 */
		private long fringedbBytes = 10*1024*1024; // 10 MiB

		/**
		 * True to print progress info to the console
		 */
		private boolean showProgress = true;

		/**
		 * Amount the threshold should be increased at each step of the sweep, in kcal/mol.
		 */
		// TODO: don't do linear steps?
		private double sweepIncrement = 1.0;

		/**
		 * How many full conformation minimizations could your hardware platform concievably do in
		 * the amount of time you're willing to wait?
		 *
		 * Used with negligableFreeEnergy to determine upper bounds for sweep thresholds.
		 *
		 * Even conformations with high energies can contibute significantly to the free energy when
		 * there are astronomical numbers of them. But if we could never possibly minimize all of them
		 * to calculate their free energy using current resources, let's just ignore those confs.
		 */
		private long maxNumMinimizations = 1000000000; // we probably can't minimize 1B confs, right?

		/**
		 * At what free energy threshold would we stop caring about the free energy contribution
		 * of a huge number of high-energy conformations?
		 *
		 * Used with maxNumMinimizations to determine upper bounds for sweep thresholds.
		 */
		private double negligableFreeEnergy = -1.0;

		public Builder(MultiStateConfSpace confSpace) {

			this.confSpace = confSpace;

			this.stateConfigs = new StateConfig[confSpace.states.size()];
		}

		public Builder configState(MultiStateConfSpace.State state, StateConfig config) {
			stateConfigs[state.index] = config;
			return this;
		}

		public Builder configEachState(Function<MultiStateConfSpace.State,StateConfig> configurator) {
			for (MultiStateConfSpace.State state : confSpace.states) {
				configState(state, configurator.apply(state));
			}
			return this;
		}

		public Builder setMathContext(MathContext val) {
			mathContext = val;
			return this;
		}

		public Builder setSeqDBFile(File val) {
			seqdbFile = val;
			return this;
		}

		public Builder setSeqDBMathContext(MathContext val) {
			seqdbMathContext = val;
			return this;
		}

		public Builder setFringeDBFile(File val) {
			fringedbFile = val;
			return this;
		}

		public Builder setFringeDBBytes(long val) {
			fringedbBytes = val;
			return this;
		}

		public Builder setFringeDBMiB(int val) {
			fringedbBytes = val*1024*1024;
			return this;
		}

		public Builder setShowProgress(boolean val) {
			showProgress = val;
			return this;
		}

		public Builder setSweepIncrement(double val) {
			sweepIncrement = val;
			return this;
		}

		public Builder setMaxNumMinimizations(long val) {
			maxNumMinimizations = val;
			return this;
		}

		public Builder setNegligableFreeEnergy(double val) {
			negligableFreeEnergy = val;
			return this;
		}

		public Sofea build() {

			// make sure all the states are configured
			List<String> unconfiguredStateNames = confSpace.states.stream()
				.filter(state -> stateConfigs[state.index] == null)
				.map(state -> state.name)
				.collect(Collectors.toList());
			if (!unconfiguredStateNames.isEmpty()) {
				throw new IllegalStateException("not all states have been configured. Please configure states: " + unconfiguredStateNames);
			}

			return new Sofea(
				confSpace,
				Arrays.asList(stateConfigs),
				mathContext,
				seqdbFile,
				seqdbMathContext,
				fringedbFile,
				fringedbBytes,
				showProgress,
				sweepIncrement,
				maxNumMinimizations,
				negligableFreeEnergy
			);
		}
	}

	/**
	 * Per-state configuration needed to run SOFEA.
	 */
	public static class StateConfig {

		/**
		 * An energy matrix of lower energy bounds on RC tuples.
		 */
		public final EnergyMatrix emat;

		/**
		 * A tool to compute energies for conformations.
		 */
		public final ConfEnergyCalculator confEcalc;

		/**
		 * File to cache conformation energies. Optional, can be null.
		 */
		public final File confDBFile;

		public StateConfig(EnergyMatrix emat, ConfEnergyCalculator confEcalc, File confDBFile) {
			this.emat = emat;
			this.confEcalc = confEcalc;
			this.confDBFile = confDBFile;
		}
	}

	public static class SeqResult {

		public final Sequence sequence;
		public final DoubleBounds lmfeFreeEnergy;
		public final DoubleBounds[] stateFreeEnergies;

		public SeqResult(Sequence sequence, DoubleBounds lmfeFreeEnergy, DoubleBounds[] stateFreeEnergies) {
			this.sequence = sequence;
			this.lmfeFreeEnergy = lmfeFreeEnergy;
			this.stateFreeEnergies = stateFreeEnergies;
		}
	}

	/** decides if computation should continue or not, and which nodes we should process */
	public static interface Criterion {

		enum Filter {
			Process,
			Requeue
		}

		enum Satisfied {
			KeepIterating,
			Terminate
		}

		default Filter filterNode(MultiStateConfSpace.State state, int[] conf, BoltzmannCalculator bcalc) {
			// accept all by default
			return Filter.Process;
		}

		Satisfied isSatisfied(SeqDB seqdb, FringeDB fringedb, long sweepCount, BoltzmannCalculator bcalc);
	}


	public final MultiStateConfSpace confSpace;
	public final List<StateConfig> stateConfigs;
	public final MathContext mathContext;
	public final File seqdbFile;
	public final MathContext seqdbMathContext;
	public final File fringedbFile;
	public final long fringedbBytes;
	public final boolean showProgress;
	public final double sweepIncrement;
	public final long maxNumMinimizations;
	public final double negligableFreeEnergy;

	public final BoltzmannCalculator bcalc;
	public final double gThresholdUpper;
	public final BigDecimal zPruneThreshold;

	private final List<StateInfo> stateInfos;

	private Sofea(MultiStateConfSpace confSpace, List<StateConfig> stateConfigs, MathContext mathContext, File seqdbFile, MathContext seqdbMathContext, File fringedbFile, long fringedbBytes, boolean showProgress, double sweepIncrement, long maxNumMinimizations, double negligableFreeEnergy) {

		this.confSpace = confSpace;
		this.stateConfigs = stateConfigs;
		this.mathContext = mathContext;
		this.seqdbFile = seqdbFile;
		this.seqdbMathContext = seqdbMathContext;
		this.fringedbFile = fringedbFile;
		this.fringedbBytes = fringedbBytes;
		this.showProgress = showProgress;
		this.sweepIncrement = sweepIncrement;
		this.maxNumMinimizations = maxNumMinimizations;
		this.negligableFreeEnergy = negligableFreeEnergy;

		bcalc = new BoltzmannCalculator(mathContext);

		gThresholdUpper = bcalc.freeEnergyPrecise(bigMath().set(bcalc.calcPrecise(negligableFreeEnergy)).div(maxNumMinimizations).get());
		zPruneThreshold = bcalc.calcPrecise(gThresholdUpper);

		// init the state info
		stateInfos = confSpace.states.stream()
			.map(state -> new StateInfo(state))
			.collect(Collectors.toList());
	}

	public StateConfig getConfig(MultiStateConfSpace.State state) {
		return stateConfigs.get(state.index);
	}

	protected StateInfo getStateInfo(MultiStateConfSpace.State state) {
		return stateInfos.get(state.index);
	}

	protected BigMath bigMath() {
		return new BigMath(mathContext);
	}

	public SeqDB openSeqDB() {
		return new SeqDB(confSpace, seqdbMathContext, seqdbFile);
	}

	public FringeDB openFringeDB() {
		if (fringedbFile.exists()) {
			return FringeDB.open(confSpace, fringedbFile);
		} else {
			return FringeDB.create(confSpace, fringedbFile, fringedbBytes, mathContext);
		}
	}

	/**
	 * Start a new design using the conf space,
	 * or fail if previous results already exist.
	 */
	public void init() {
		init(false);
	}

	/**
	 * Start a new design using the conf space, or fail unless overwrite=true.
	 * If overwrite=true, the previous results are destroyed.
	 */
	public void init(boolean overwrite) {

		// don't overwrite an existing results unless explicitly asked
		if (!overwrite && (seqdbFile.exists() || fringedbFile.exists())) {
			throw new Error("Database files already exist, will not destroy existing results at:"
				+ "\n\t" + seqdbFile.getAbsolutePath()
				+ "\n\t" + fringedbFile.getAbsolutePath()
			);
		}

		// clear old results
		seqdbFile.delete();
		fringedbFile.delete();

		try (SeqDB seqdb = openSeqDB()) {
			try (FringeDB fringedb = openFringeDB()) {

				// process the root node for each state
				FringeDB.Transaction fringetx = fringedb.transaction();
				SeqDB.Transaction seqtx = seqdb.transaction();
				for (MultiStateConfSpace.State state : confSpace.states) {
					StateInfo stateInfo = stateInfos.get(state.index);

					// get a multi-sequence Z bound on the root node
					ConfIndex index = stateInfo.makeConfIndex();
					BigDecimal zSumUpper = stateInfo.calcZSumUpper(index, stateInfo.rcs);

					// init the fringe with the root node
					fringetx.writeRootNode(state, zSumUpper);
					seqtx.addZSumUpper(state, state.confSpace.makeUnassignedSequence(), zSumUpper);
				}
				fringetx.commit();
				fringedb.finishSweep();
				seqtx.commit();
			}
		}
	}

	private static class Node {

		final int[] conf;
		final BigDecimal zSumUpper;

		Node(int[] conf, BigDecimal zSumUpper) {
			this.conf = conf;
			this.zSumUpper = zSumUpper;
		}
	}

	private static class ZPath {

		final int[] conf;
		final BigDecimal zPath;

		ZPath(int[] conf, BigDecimal zPath) {
			this.conf = conf;
			this.zPath = zPath;
		}
	}

	private class NodeTransaction {

		final MultiStateConfSpace.State state;
		final int[] conf;
		final BigDecimal zSumUpper;

		final ConfIndex index;
		final List<Node> replacementNodes = new ArrayList<>();
		final List<ZPath> zPaths = new ArrayList<>();

		NodeTransaction(MultiStateConfSpace.State state, int[] conf, BigDecimal zSumUpper) {

			this.state = state;
			this.conf = conf;
			this.zSumUpper = zSumUpper;

			index = new ConfIndex(state.confSpace.positions.size());
			Conf.index(conf, index);
		}

		void addReplacementNode(ConfIndex index, BigDecimal zSumUpper) {
			replacementNodes.add(new Node(Conf.make(index), zSumUpper));
		}

		int numReplacementNodes() {
			return replacementNodes.size();
		}

		void addZPath(ConfIndex index, BigDecimal zPath) {
			zPaths.add(new ZPath(Conf.make(index), zPath));
		}

		int  numZPaths() {
			return zPaths.size();
		}

		boolean hasRoomToReplace(FringeDB.Transaction fringetx, int otherNodesInFlight) {
			return fringetx.dbHasRoomFor(replacementNodes.size() + otherNodesInFlight);
		}

		boolean replace(FringeDB.Transaction fringetx, SeqDB.Transaction seqtx) {

			// flush transactions if needed
			boolean flush = !fringetx.txHasRoomFor(replacementNodes.size());
			if (flush) {
				flushTransactions(fringetx, seqtx);
			}

			StateInfo stateInfo = stateInfos.get(state.index);

			// subtract the node zSumUpper
			seqtx.subZSumUpper(state, stateInfo.makeSeq(conf), zSumUpper);

			// update fringedb and seqdb with the replacement nodes
			for (Node replacementNode : replacementNodes) {
				fringetx.writeReplacementNode(state, replacementNode.conf, replacementNode.zSumUpper);
				seqtx.addZSumUpper(state, stateInfo.makeSeq(replacementNode.conf), replacementNode.zSumUpper);
			}

			// add the z values
			for (ZPath zPath : zPaths) {
				seqtx.addZPath(state, stateInfo.makeSeq(zPath.conf), zPath.zPath);
			}

			return flush;
		}

		boolean requeue(FringeDB.Transaction fringetx, SeqDB.Transaction seqtx) {

			// flush transactions if needed
			boolean flush = !fringetx.txHasRoomFor(1);
			if (flush) {
				flushTransactions(fringetx, seqtx);
			}

			// ignore all of the seqdb changes

			// move the node to the end of the fringedb queue
			fringetx.writeReplacementNode(state, conf, zSumUpper);

			return flush;
		}

		boolean flushTransactions(FringeDB.Transaction fringetx, SeqDB.Transaction seqtx) {

			// make sure we didn't overflow the buffer entirely
			if (replacementNodes.size() > fringetx.maxWriteBufferNodes()) {
				throw new IllegalStateException(String.format("FringeDB write buffer is too small. Holds %d nodes, but need %d nodes",
					fringetx.maxWriteBufferNodes(), replacementNodes.size()
				));
			}

			// commit both transactions at the same time
			fringetx.commit();
			seqtx.commit();

			return true;
		}
	}

	private class ConfTables implements AutoCloseable {

		StateInfo.Confs[] confsByState;

		ConfTables() {
			confsByState = new StateInfo.Confs[confSpace.states.size()];
			for (MultiStateConfSpace.State state : confSpace.states) {
				confsByState[state.index] = stateInfos.get(state.index).new Confs();
			}
		}

		@Override
		public void close() {
			for (StateInfo.Confs confs : confsByState) {
				confs.close();
			}
		}

		public ConfDB.ConfTable get(MultiStateConfSpace.State state) {
			return confsByState[state.index].table;
		}
	}

	private static class StateStats {
		long read = 0;
		long expanded = 0;
		long requeuedByThreshold = 0;
		long requeuedByFilter = 0;
		long requeuedForSpace = 0;
		long added = 0;
		long minimizations = 0;
	}

	/**
	 * Keep sweeping the threshold over the parameter space until the criterion is satisfied.
	 */
	public void refine(Criterion criterion) {
		try (SeqDB seqdb = openSeqDB()) {
		try (FringeDB fringedb = openFringeDB()) {
		try (ConfTables confTables = new ConfTables()) {

			// use the energy calculator to provide the paralleism
			TaskExecutor tasks = stateConfigs.get(0).confEcalc.tasks;

			// calculate lower bounds on gThreshold
			double[] gThresholdLower = confSpace.states.stream()
				.mapToDouble(state -> {
					Sofea.StateInfo stateInfo = getStateInfo(state);
					BigDecimal zSumUpper = stateInfo.calcZSumUpper(stateInfo.makeConfIndex(), stateInfo.rcs);
					return bcalc.freeEnergyPrecise(zSumUpper);
				})
				.toArray();

			// init the gThresholds based on the fringe set
			Double[] gThreshold = confSpace.states.stream()
				.map(state -> {
					BigDecimal zSumMax = fringedb.getZSumMax(state);
					if (zSumMax == null) {
						return null;
					} else {
						return bcalc.freeEnergyPrecise(zSumMax);
					}
				})
				.toArray(size -> new Double[size]);

			for (long step=1; ; step++) {

				// check the termination criterion
				if (criterion != null && criterion.isSatisfied(seqdb, fringedb, step, bcalc) == Criterion.Satisfied.Terminate) {
					log("SOFEA finished, criterion satisfied");
					break;
				}

				// stop if we ran out of fringe nodes
				if (fringedb.isEmpty()) {
					log("SOFEA finished, explored every node");
					break;
				}

				// increment gThreshold before each pass over the fringe set
				for (MultiStateConfSpace.State state : confSpace.states) {
					if (gThreshold[state.index] != null) {
						if (gThreshold[state.index] > gThresholdUpper) {
							gThreshold[state.index] = null;
						} else {
							gThreshold[state.index] += sweepIncrement;
						}
					}
				}

				// calc the corresponding zThreshold
				BigDecimal[] zThreshold = confSpace.states.stream()
					.map(state -> {
						Double g = gThreshold[state.index];
						if (g == null) {
							return null;
						} else {
							return bcalc.calcPrecise(g);
						}
					})
					.toArray(size -> new BigDecimal[size]);

				// show progress if needed
				if (showProgress) {
					log("step %d", step);
					for (MultiStateConfSpace.State state : confSpace.states) {
						log("\t%10s  gThreshold = %9.3f  in  [%9.3f,%9.3f]",
							state.name,
							gThreshold[state.index],
							gThresholdLower[state.index],
							gThresholdUpper
						);
					}
				}

				// init sweep stats
				StateStats[] stats = new StateStats[confSpace.states.size()];
				for (MultiStateConfSpace.State state : confSpace.states) {
					stats[state.index] = new StateStats();
				}

				// keep track of how many nodes are in outstanding tasks, and hence unknown to FringeDB's size counters
				// NOTE: use a size-one array instead of a plain var, since Java's compiler is kinda dumb about lambdas
				int[] nodesInFlight = { 0 };

				// start (or resume) the iteration over the fringe set
				FringeDB.Transaction fringetx = fringedb.transaction();
				SeqDB.Transaction seqtx = seqdb.transaction();
				Progress progress = new Progress(fringetx.numNodesToRead());
				progress.setReportMemory(true);
				while (true) {

					// read the next node and make a transaction for it
					final NodeTransaction nodetx;
					synchronized (Sofea.this) { // don't race the listener thread
						if (!fringetx.hasNodesToRead()) {
							break;
						}
						fringetx.readNode();
						nodesInFlight[0]++;
						nodetx = new NodeTransaction(
							fringetx.state(),
							fringetx.conf(),
							fringetx.zSumUpper()
						);
					}
					stats[nodetx.state.index].read++;

					// check the node filter in the criterion
					if (criterion.filterNode(nodetx.state, nodetx.conf, bcalc) == Criterion.Filter.Requeue) {
						synchronized (Sofea.this) { // don't race the listener thread
							stats[nodetx.state.index].requeuedByFilter++;
							nodesInFlight[0]--;
							nodetx.requeue(fringetx, seqtx);
							if (showProgress) {
								progress.incrementProgress();
							}
						}
						continue;
					}

					// process nodes with tasks (possibly in parallel)
					tasks.submit(
						() -> design(
							nodetx,
							zThreshold[nodetx.state.index],
							nodetx.index,
							nodetx.zSumUpper,
							confTables.get(nodetx.state)
						),
						(wasExpanded) -> {

							stats[nodetx.state.index].minimizations += nodetx.numZPaths();
							if (wasExpanded) {

								synchronized (Sofea.this) { // don't race the main thread
									nodesInFlight[0]--;
									if (nodetx.hasRoomToReplace(fringetx, nodesInFlight[0])) {
										stats[nodetx.state.index].added += nodetx.numReplacementNodes();
										stats[nodetx.state.index].expanded++;
										nodetx.replace(fringetx, seqtx);
									} else {
										stats[nodetx.state.index].requeuedForSpace++;
										nodetx.requeue(fringetx, seqtx);
									}
									if (showProgress) {
										progress.incrementProgress();
									}
								}

							} else {

								stats[nodetx.state.index].requeuedByThreshold++;
								synchronized (Sofea.this) { // don't race the main thread
									nodesInFlight[0]--;
									nodetx.requeue(fringetx, seqtx);
									if (showProgress) {
										progress.incrementProgress();
									}
								}
							}
						}
					);
				}
				tasks.waitForFinish();
				fringetx.commit();
				seqtx.commit();
				fringedb.finishSweep();

				// show stats if needed
				if (showProgress) {
					log("\tfringe size: %d/%d (%.1f%%) nodes",
						fringedb.getNumNodes(), fringedb.getCapacity(),
						100.0f*fringedb.getNumNodes()/fringedb.getCapacity()
					);
					for (MultiStateConfSpace.State state : confSpace.states) {
						log("\t%10s  read=%6d [ expanded=%6d  requeuedByThreshold=%6d  requeuedByFilter=%6d  requeuedForSpace=%6d ] added=%6d  minimizations=%6d",
							state.name,
							stats[state.index].read,
							stats[state.index].expanded,
							stats[state.index].requeuedByThreshold,
							stats[state.index].requeuedByFilter,
							stats[state.index].requeuedForSpace,
							stats[state.index].added,
							stats[state.index].minimizations
						);
					}
					log("\tNon-garbage heap memory usage: %s", JvmMem.getOldPool());
				}
			}
		}}}
	}

	private boolean design(NodeTransaction nodetx, BigDecimal zSumThreshold, ConfIndex index, BigDecimal zSumUpper, ConfDB.ConfTable confTable) {

		// forget any subtree if we've hit the end of the sweep interval, or that's below the pruning threshold
		if (zSumThreshold == null || MathTools.isLessThan(zSumUpper, zPruneThreshold)) {
			return false;
		}

		// if zSumUpper is too small, add the node to the fringe set
		if (MathTools.isLessThan(zSumUpper, zSumThreshold)) {
			nodetx.addReplacementNode(index, zSumUpper);
			return false;
		}

		StateInfo stateInfo = stateInfos.get(nodetx.state.index);

		// is this a leaf node?
		if (index.isFullyDefined()) {

			// yup, add the zPath
			nodetx.addZPath(index, stateInfo.calcZPath(index, confTable));
			return true;
		}

		// not a leaf node

		// no more easy bound improvement possible, so recuse
		int pos = stateInfo.posPermutation[index.numDefined];
		for (int rc : stateInfo.rcs.get(pos)) {
			index.assignInPlace(pos, rc);
			design(
				nodetx,
				zSumThreshold,
				index,
				stateInfo.calcZSumUpper(index, stateInfo.rcs),
				confTable
			);
			index.unassignInPlace(pos);
		}

		return true;
	}

	/** WARNING: naive brute force method, for testing small trees only */
	protected BigDecimal calcZSum(Sequence seq, MultiStateConfSpace.State state) {
		StateInfo stateInfo = stateInfos.get(state.index);
		try (StateInfo.Confs confs = stateInfo.new Confs()) {
			RCs rcs = seq.makeRCs(state.confSpace);
			return stateInfo.calcZSum(stateInfo.makeConfIndex(), rcs, confs.table);
		}
	}

	protected class StateInfo {

		class Confs implements AutoCloseable {

			final ConfDB confdb;
			final ConfDB.ConfTable table;

			Confs() {
				StateConfig config = stateConfigs.get(state.index);
				if (config.confDBFile != null) {
					confdb = new ConfDB(state.confSpace, config.confDBFile);
					table = confdb.new ConfTable("sofea");
				} else {
					confdb = null;
					table = null;
				}
			}

			@Override
			public void close() {
				if (confdb != null) {
					confdb.close();
				}
			}
		}

		final MultiStateConfSpace.State state;

		final ConfEnergyCalculator confEcalc;
		final ZMatrix zmat;
		final RCs rcs;
		final int[] posPermutation;

		private final int[][] rtsByRcByPos;
		private final int[] numRtsByPos;

		StateInfo(MultiStateConfSpace.State state) {

			this.state = state;
			StateConfig config = stateConfigs.get(state.index);
			this.confEcalc = config.confEcalc;
			this.zmat = new ZMatrix(state.confSpace);
			this.zmat.set(config.emat, bcalc);
			this.rcs = new RCs(state.confSpace);

			// sort positions so multi-sequence layers are first
			posPermutation = state.confSpace.positions.stream()
				.sorted((a, b) -> {

					// prefer mutable positions first
					if (a.hasMutations() && !b.hasMutations()) {
						return -1;
					} else if (!a.hasMutations() && b.hasMutations()) {
						return +1;
					}

					/* TODO: what heuristic works well here?
					// then, sort by ordering heuristic
					// negate to sort in (weakly) descending order
					int order = -calcOrderHeuristic(a).compareTo(calcOrderHeuristic(b));
					if (order != 0) {
						return order;
					}
					*/

					// otherwise, sort by index
					return a.index - b.index;
				})
				.mapToInt(pos -> pos.index)
				.toArray();

			// calculate all the RTs by RC and pos
			rtsByRcByPos = new int[state.confSpace.positions.size()][];
			numRtsByPos = new int[state.confSpace.positions.size()];
			for (SimpleConfSpace.Position confPos : state.confSpace.positions) {
				rtsByRcByPos[confPos.index] = new int[confPos.resConfs.size()];
				SeqSpace.Position seqPos = confSpace.seqSpace.getPosition(confPos.resNum);
				if (seqPos != null) {
					numRtsByPos[confPos.index] = seqPos.resTypes.size();
					for (SimpleConfSpace.ResidueConf rc : confPos.resConfs) {
						SeqSpace.ResType rt = seqPos.getResTypeOrThrow(rc.template.name);
						rtsByRcByPos[confPos.index][rc.index] = rt.index;
					}
				} else {
					numRtsByPos[confPos.index] = 0;
					Arrays.fill(rtsByRcByPos[confPos.index], Sequence.Unassigned);
				}
			}
		}

		ConfIndex makeConfIndex() {
			ConfIndex index = new ConfIndex(state.confSpace.positions.size());
			index.updateUndefined();
			return index;
		}

		Sequence makeSeq(int[] conf) {
			return confSpace.seqSpace.makeSequence(state.confSpace, conf);
		}

		/** WARNING: naive brute force method, for testing small trees only */
		Map<Sequence,BigInteger> calcNumLeavesBySequence(ConfIndex index) {

			Map<Sequence,BigInteger> out = new HashMap<>();

			// NOTE: java has a hard time with recursive lambdas,
			// so use an array to work around the compiler's limitations
			Runnable[] f = { null };
			f[0] = () -> {

				// if this is a leaf node, increment the sequence counter
				if (index.isFullyDefined()) {
					out.compute(makeSeq(Conf.make(index)), (seq,  old) -> {
						if (old == null) {
							return BigInteger.ONE;
						} else {
							return old.add(BigInteger.ONE);
						}
					});
					return;
				}

				// otherwise, recurse
				int pos = posPermutation[index.numDefined];
				for (int rc : rcs.get(pos)) {
					index.assignInPlace(pos, rc);
					f[0].run();
					index.unassignInPlace(pos);
				}
			};
			f[0].run();

			return out;
		}

		BigInteger calcNumLeavesUpperBySequence(ConfIndex index, RCs rcs) {

			BigInteger count = BigInteger.ONE;

			for (int i=0; i<index.numUndefined; i++) {
				int pos = index.undefinedPos[i];

				int maxCount = 0;

				int numRts = numRtsByPos[pos];
				if (numRts > 0) {

					// count the RCs by RT
					int[] counts = new int[numRts];
					for (int rc : rcs.get(pos)) {
						int rtCount = ++counts[rtsByRcByPos[pos][rc]];
						maxCount = Math.max(maxCount, rtCount);
					}

				} else {

					// all RCs are the same RT
					maxCount = rcs.getNum(pos);
				}

				count = count.multiply(BigInteger.valueOf(maxCount));
			}

			return count;
		}

		BigDecimal calcZSumUpper(ConfIndex index, RCs rcs) {
			return bigMath()
				.set(calcZPathHeadUpper(index))
				.mult(calcZPathTailUpper(index, rcs))
				.mult(calcNumLeavesUpperBySequence(index, rcs))
				.get();
		}

		BigDecimal calcZPathUpper(ConfIndex index, RCs rcs) {
			return bigMath()
				.set(calcZPathHeadUpper(index))
				.mult(calcZPathTailUpper(index, rcs))
				.get();
		}

		BigDecimal calcZPathHeadUpper(ConfIndex index) {

			BigMath z = bigMath().set(1.0);

			// multiply all the singles and pairs
			for (int i1=0; i1<index.numDefined; i1++) {
				int pos1 = index.definedPos[i1];
				int rc1 = index.definedRCs[i1];

				z.mult(zmat.getOneBody(pos1, rc1));

				for (int i2=0; i2<i1; i2++) {
					int pos2 = index.definedPos[i2];
					int rc2 = index.definedRCs[i2];

					z.mult(zmat.getPairwise(pos1, rc1, pos2, rc2));
				}
			}

			// multiply the higher order corrections if needed
			zmat.forEachHigherOrderTupleIn(index, (tuple, tupleZ) -> {
				z.mult(tupleZ);
			});

			return z.get();
		}

		BigDecimal calcZPathTailUpper(ConfIndex index, RCs rcs) {

			// TODO: use higher-order corrections?

			// this is the usual A* heuristic

			MathTools.Optimizer opt = MathTools.Optimizer.Maximize;
			BigMath z = bigMath().set(1.0);

			// for each undefined position
			for (int i=0; i<index.numUndefined; i++) {
				int pos1 = index.undefinedPos[i];

				// optimize over possible assignments to pos1
				BigDecimal zpos1 = opt.initBigDecimal();
				for (int rc1 : rcs.get(pos1)) {

					BigMath zrc1 = bigMath().set(zmat.getOneBody(pos1, rc1));

					// interactions with defined residues
					for (int j=0; j<index.numDefined; j++) {
						int pos2 = index.definedPos[j];
						int rc2 = index.definedRCs[j];

						zrc1.mult(zmat.getPairwise(pos1, rc1, pos2, rc2));
					}

					// interactions with undefined residues
					for (int j=0; j<i; j++) {
						int pos2 = index.undefinedPos[j];

						// optimize over possible assignments to pos2
						// TODO: make this a look-up table?
						BigDecimal optzrc2 = opt.initBigDecimal();
						for (int rc2 : rcs.get(pos2)) {

							// pair with pos2
							BigDecimal zrc2 = zmat.getPairwise(pos1, rc1, pos2, rc2);

							optzrc2 = opt.opt(optzrc2, zrc2);
						}

						zrc1.mult(optzrc2);
					}

					zpos1 = opt.opt(zpos1, zrc1.get());
				}

				assert (MathTools.isFinite(zpos1));
				z.mult(zpos1);
			}

			return z.get();
		}

		BigDecimal calcZPath(ConfIndex index, ConfDB.ConfTable confTable) {

			if (!index.isFullyDefined()) {
				throw new IllegalArgumentException("not a full conf");
			}

			double e = confEcalc.calcEnergy(new RCTuple(index), confTable);
			return bcalc.calcPrecise(e);
		}

		/** a reminder to myself this this is a bad idea */
		@Deprecated
		BigDecimal calcZPathHead(ConfIndex index, ConfDB.ConfTable confTable) {
			throw new Error("stop trying to do this! There's no such thing!");
		}

		/** WARNING: naive brute force method, for testing small trees only */
		BigDecimal calcZSum(ConfIndex index, RCs rcs, ConfDB.ConfTable confTable) {

			// base case, compute the conf energy
			if (index.isFullyDefined()) {
				double e = confEcalc.calcEnergy(new RCTuple(index), confTable);
				return bcalc.calcPrecise(e);
			}

			// otherwise, recurse
			BigMath m = bigMath().set(0.0);
			int pos = posPermutation[index.numDefined];
			for (int rc : rcs.get(pos)) {
				index.assignInPlace(pos, rc);
				m.add(calcZSum(index, rcs, confTable));
				index.unassignInPlace(pos);
			}
			return m.get();
		}

		/** WARNING: naive brute force method, for testing small trees only */
		BigDecimalBounds calcZPathBoundsExact(ConfIndex index, RCs rcs, ConfDB.ConfTable confTable) {

			// base case, compute the conf energy
			if (index.isFullyDefined()) {
				return new BigDecimalBounds(calcZPath(index, confTable));
			}

			// otherwise, recurse
			BigMath mlo = bigMath();
			BigMath mhi = bigMath();
			int pos = posPermutation[index.numDefined];
			for (int rc : rcs.get(pos)) {
				index.assignInPlace(pos, rc);
				BigDecimalBounds sub = calcZPathBoundsExact(index, rcs, confTable);
				mlo.minOrSet(sub.lower);
				mhi.maxOrSet(sub.upper);
				index.unassignInPlace(pos);
			}
			return new BigDecimalBounds(mlo.get(), mhi.get());
		}
	}
}
