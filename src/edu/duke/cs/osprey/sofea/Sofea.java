package edu.duke.cs.osprey.sofea;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.*;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import edu.duke.cs.osprey.tools.MathTools.BigIntegerBounds;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;


/**
 * TODO: doc me!
 */
public class Sofea {

	public static class Builder {

		public final MultiStateConfSpace confSpace;

		private StateConfig[] stateConfigs;
		private MathContext mathContext = new MathContext(32, RoundingMode.HALF_UP);
		private File seqdbFile = new File("seq.db");
		private MathContext seqdbMathContext = new MathContext(128, RoundingMode.HALF_UP);
		private File fringedbFile = new File("fringe.db");
		private long fringedbBytes = 10*1024*1024; // 10 MiB
		private boolean showProgress = true;
		private Parallelism parallelism = Parallelism.makeCpu(1);
		private double sweepDivisor = Math.pow(Math.E, 4.0);
		private double zPruneThreshold = 1e-5; // ln(1 + 1e-5) < 0.0000

		// NOTE: don't need much precision for most math, but need lots of precision for seqdb math

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

		public Builder setParallelism(Parallelism val) {
			parallelism = val;
			return this;
		}

		public Builder setSweepDivisor(double val) {
			sweepDivisor = val;
			return this;
		}

		public Builder setZPruneThreshold(double val) {
			zPruneThreshold = val;
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
				parallelism,
				sweepDivisor,
				zPruneThreshold
			);
		}
	}

	public static class StateConfig {

		public final EnergyMatrix ematLower;
		public final EnergyMatrix ematUpper;
		public final ConfEnergyCalculator confEcalc;
		public final File confDBFile;

		public StateConfig(EnergyMatrix ematLower, EnergyMatrix ematUpper, ConfEnergyCalculator confEcalc, File confDBFile) {
			this.ematLower = ematLower;
			this.ematUpper = ematUpper;
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
	public final Parallelism parallelism;
	public final double sweepDivisor;
	public final BigDecimal zPruneThreshold;

	public final BoltzmannCalculator bcalc;

	private final List<StateInfo> stateInfos;

	private Sofea(MultiStateConfSpace confSpace, List<StateConfig> stateConfigs, MathContext mathContext, File seqdbFile, MathContext seqdbMathContext, File fringedbFile, long fringedbBytes, boolean showProgress, Parallelism parallelism, double sweepDivisor, double zPruneThreshold) {

		this.confSpace = confSpace;
		this.stateConfigs = stateConfigs;
		this.mathContext = mathContext;
		this.seqdbFile = seqdbFile;
		this.seqdbMathContext = seqdbMathContext;
		this.fringedbFile = fringedbFile;
		this.fringedbBytes = fringedbBytes;
		this.showProgress = showProgress;
		this.parallelism = parallelism;
		this.sweepDivisor = sweepDivisor;
		this.zPruneThreshold = MathTools.biggen(zPruneThreshold);

		bcalc = new BoltzmannCalculator(mathContext);

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
	 * start a new design using the conf space,
	 * or fail if previous results already exist
	 */
	public void init() {
		init(false);
	}

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
					BigDecimalBounds zSumBounds = stateInfo.calcZSumBounds(index, stateInfo.rcs);

					// init the fringe with the root node
					fringetx.writeRootNode(state, zSumBounds, new BigDecimalBounds(BigDecimal.ONE));
					seqtx.addZSumBounds(state, state.confSpace.makeUnassignedSequence(), zSumBounds);
				}
				fringetx.commit();
				fringedb.finishSweep();
				seqtx.commit();
			}
		}
	}

	private static class Node {

		final int[] conf;
		final BigDecimalBounds zSumBounds;
		final BigDecimalBounds zPathHeadBounds;

		Node(int[] conf, BigDecimalBounds zSumBounds, BigDecimalBounds zPathHeadBounds) {
			this.conf = conf;
			this.zSumBounds = zSumBounds;
			this.zPathHeadBounds = zPathHeadBounds;
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
		final BigDecimalBounds zSumBounds;
		final BigDecimalBounds zPathHeadBounds;

		final ConfIndex index;
		final List<Node> replacementNodes = new ArrayList<>();
		final List<ZPath> zPaths = new ArrayList<>();

		NodeTransaction(MultiStateConfSpace.State state, int[] conf, BigDecimalBounds zSumBounds, BigDecimalBounds zPathHeadBounds) {

			this.state = state;
			this.conf = conf;
			this.zSumBounds = zSumBounds;
			this.zPathHeadBounds = zPathHeadBounds;

			index = new ConfIndex(state.confSpace.positions.size());
			Conf.index(conf, index);
		}

		void addReplacementNode(ConfIndex index, BigDecimalBounds zSumBounds, BigDecimalBounds zPathHeadBounds) {
			replacementNodes.add(new Node(Conf.make(index), zSumBounds, zPathHeadBounds));
		}

		int numReplacementNodes() {
			return replacementNodes.size();
		}

		void addZPath(ConfIndex index, BigDecimal zPath) {
			zPaths.add(new ZPath(Conf.make(index), zPath));
		}

		boolean hasRoomToReplace(FringeDB.Transaction fringetx, int otherNodesInFlight) {
			return fringetx.dbHasRoomFor(replacementNodes.size() + otherNodesInFlight);
		}

		boolean replace(FringeDB.Transaction fringetx, SeqDB.Transaction seqtx) {

			boolean flushed = flushTransactionsIfNeeded(fringetx, seqtx);

			StateInfo stateInfo = stateInfos.get(state.index);

			// subtract the node zSumBounds
			seqtx.subZSumBounds(state, stateInfo.makeSeq(conf), zSumBounds);

			// update fringedb and seqdb with the replacement nodes
			for (Node replacementNode : replacementNodes) {
				fringetx.writeReplacementNode(state, replacementNode.conf, replacementNode.zSumBounds, replacementNode.zPathHeadBounds);
				seqtx.addZSumBounds(state, stateInfo.makeSeq(replacementNode.conf), replacementNode.zSumBounds);
			}

			// add the z values
			for (ZPath zPath : zPaths) {
				seqtx.addZPath(state, stateInfo.makeSeq(zPath.conf), zPath.zPath);
			}

			return flushed;
		}

		boolean requeue(FringeDB.Transaction fringetx, SeqDB.Transaction seqtx) {

			boolean flushed = flushTransactionsIfNeeded(fringetx, seqtx);

			// ignore all of the seqdb changes

			// move the node to the end of the fringedb queue
			fringetx.writeReplacementNode(state, conf, zSumBounds, zPathHeadBounds);

			return flushed;
		}

		boolean flushTransactionsIfNeeded(FringeDB.Transaction fringetx, SeqDB.Transaction seqtx) {

			// skip if the fringe db transaction isn't full
			if (fringetx.txHasRoomFor(replacementNodes.size())) {
				return false;
			}

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
		long replaced = 0;
		long added = 0;
		long requeued = 0;
	}

	/**
	 * keep sweeping over the threshold space until the criterion is satisfied
	 */
	public void refine(Criterion criterion) {
		try (SeqDB seqdb = openSeqDB()) {
		try (FringeDB fringedb = openFringeDB()) {
		try (ConfTables confTables = new ConfTables()) {
		try (TaskExecutor tasks = parallelism.makeTaskExecutor()) {

			// init the zSumWidthThresholds
			BigDecimal[] zSumWidthThreshold = confSpace.states.stream()
				.map(state -> fringedb.getZSumMax(state))
				.toArray(size -> new BigDecimal[size]);

			long sweepCount = 0;
			while (true) {

				// check the termination criterion
				if (criterion != null && criterion.isSatisfied(seqdb, fringedb, sweepCount, bcalc) == Criterion.Satisfied.Terminate) {
					log("SOFEA finished, criterion satisfied");
					break;
				}

				// stop if we ran out of fringe nodes
				if (fringedb.isEmpty()) {
					log("SOFEA finished, explored every node");
					break;
				}

				// reduce zSumWidthThreshold before each sweep
				for (MultiStateConfSpace.State state : confSpace.states) {
					BigDecimal zSumMax = fringedb.getZSumMax(state);
					if (zSumMax == null) {
						zSumWidthThreshold[state.index] = null;
					} else {
						zSumWidthThreshold[state.index] = bigMath()
							.set(zSumWidthThreshold[state.index])
							.min(zSumMax)
							.div(sweepDivisor)
							.get();
					}
				}

				// show progress if needed
				sweepCount++;
				if (showProgress) {
					log("step %d", sweepCount);
					for (MultiStateConfSpace.State state : confSpace.states) {
						if (zSumWidthThreshold[state.index] == null) {
							log("\t%10s  finished", state.name);
						} else {
							log("\t%10s  ln1p(Z sum width threshold) = %s/%s",
								state.name,
								Log.formatBigLn(zSumWidthThreshold[state.index]),
								Log.formatBigLn(fringedb.getZSumMax(state))
							);
						}
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

				// start (or resume) the sweep
				FringeDB.Transaction fringetx = fringedb.transaction();
				SeqDB.Transaction seqtx = seqdb.transaction();
				Progress progress = new Progress(fringetx.numNodesToRead());
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
							fringetx.zSumBounds(),
							fringetx.zPathHeadBounds()
						);
					}
					stats[nodetx.state.index].read++;

					// check the node filter in the criterion
					if (criterion.filterNode(fringetx.state(), fringetx.conf(), bcalc) == Criterion.Filter.Requeue) {
						stats[nodetx.state.index].requeued++;
						nodetx.requeue(fringetx, seqtx);
						continue;
					}

					// process nodes with tasks (possibly in parallel)
					tasks.submit(
						() -> design(
							nodetx,
							zSumWidthThreshold[nodetx.state.index],
							nodetx.index,
							nodetx.zSumBounds,
							nodetx.zPathHeadBounds,
							confTables.get(nodetx.state)
						),
						(wasExpanded) -> {

							if (wasExpanded) {
								stats[nodetx.state.index].expanded++;
							}

							synchronized (Sofea.this) { // don't race the main thread
								nodesInFlight[0]--;
								if (nodetx.hasRoomToReplace(fringetx, nodesInFlight[0])) {
									stats[nodetx.state.index].added += nodetx.numReplacementNodes();
									stats[nodetx.state.index].replaced++;
									nodetx.replace(fringetx, seqtx);
								} else {
									stats[nodetx.state.index].requeued++;
									nodetx.requeue(fringetx, seqtx);
								}
							}

							if (showProgress) {
								progress.incrementProgress();
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
						log("\t%10s  read=%6d  expanded=%6d  replaced=%6d  added=%6d requeued=%6d",
							state.name,
							stats[state.index].read,
							stats[state.index].expanded,
							stats[state.index].replaced,
							stats[state.index].added,
							stats[state.index].requeued
						);
					}
				}
			}
		}}}}
	}

	private boolean design(NodeTransaction nodetx, BigDecimal zSumWidthThreshold, ConfIndex index, BigDecimalBounds zSumBounds, BigDecimalBounds zPathHeadBounds, ConfDB.ConfTable confTable) {

		// forget any subtree that's below the pruning threshold
		/* NOTE:
			Unfortunately, a hard pruning threshold is necessary in practice.
			Clashing confs can have extremely low zPath values (like 1e-30,000,000).
			When FringeDB runs out of space, we'd normally have to wait for zSumWidthThreshold
			to drop below the zPath of all leaves in a subtree before we can remove that whole
			subtree from FringeDB. But if the lowest zPath is supremely low, zSumWidthThreshold
			will never drop that low even after millions of steps in the sweep. So we need to
			just give up on these subtrees entirely and hope astronomical numbers of those confs
			don't add up to a significant portion of the zSum.
		 */
		if (MathTools.isLessThan(zSumBounds.upper, zPruneThreshold)) {
			return false;
		}

		// if zSumBounds width is too small, add the node to the fringe set
		if (MathTools.isLessThan(zSumBounds.size(mathContext), zSumWidthThreshold)) {
			nodetx.addReplacementNode(index, zSumBounds, zPathHeadBounds);
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

		/* TODO: the extra work to improve the bounds never seems to pay off, maybe there's a better way?
		// maybe we do more work to improve zSumBounds and that will reduce uncertainty without recursing?
		// more work is expensive, so only do it if we're high up the tree and hence the payoff could be large
		if (index.numUndefined >= 5) {

			// TODO: and for muli-sequence nodes?
			if (stateInfo.isSingleSequence(index)) {

				// don't overwrite the existing bounds, make a new ones
				BigDecimalBounds zPathHeadBoundsTighter = stateInfo.calcZPathHeadBoundsTighter(index, confTable);
				BigDecimalBounds zSumBoundsTighter = stateInfo.calcZSumBoundsTighter(index, stateInfo.rcs, confTable, zPathHeadBoundsTighter);

				// if the width is too small now, add it to the fringe set
				if (MathTools.isLessThan(zSumBounds.size(mathContext), zSumWidthThreshold)) {
					nodetx.addReplacementNode(index, zSumBoundsTighter, zPathHeadBoundsTighter);
					return false;
				}
			}
		}
		*/

		// no more bound improvement possible, recuse
		int pos = stateInfo.posPermutation[index.numDefined];
		for (int rc : stateInfo.rcs.get(pos)) {

			BigDecimalBounds zPathHeadBoundsRc = multBounds(
				zPathHeadBounds,
				stateInfo.calcZPathNodeBounds(index, pos, rc)
			);

			index.assignInPlace(pos, rc);

			BigDecimalBounds zSumBoundsRc = stateInfo.calcZSumBounds(index, stateInfo.rcs);
			design(
				nodetx,
				zSumWidthThreshold,
				index,
				zSumBoundsRc,
				zPathHeadBoundsRc,
				confTable
			);

			index.unassignInPlace(pos);
		}

		return true;
	}

	/** WARNING: naive brute force method, for testing small trees only */
	public BigDecimal calcZSum(Sequence seq, MultiStateConfSpace.State state) {
		StateInfo stateInfo = stateInfos.get(state.index);
		try (StateInfo.Confs confs = stateInfo.new Confs()) {
			RCs rcs = seq.makeRCs(state.confSpace);
			return stateInfo.calcZSum(stateInfo.makeConfIndex(), rcs, confs.table);
		}
	}

	protected class StateInfo {

		class Confs implements AutoCloseable {

			ConfDB confdb;
			ConfDB.ConfTable table;

			Confs() {
				StateConfig config = stateConfigs.get(state.index);
				confdb = new ConfDB(state.confSpace, config.confDBFile);
				table = confdb.new ConfTable("sofea");
			}

			@Override
			public void close() {
				confdb.close();
			}
		}

		final MultiStateConfSpace.State state;

		final ConfEnergyCalculator confEcalc;
		final ZMatrix zmatLower;
		final ZMatrix zmatUpper;
		final RCs rcs;
		final int[] posPermutation;

		private final int[][] rtsByRcByPos;
		private final int[] numRtsByPos;

		StateInfo(MultiStateConfSpace.State state) {

			this.state = state;
			StateConfig config = stateConfigs.get(state.index);
			this.confEcalc = config.confEcalc;
			this.zmatLower = new ZMatrix(state.confSpace);
			this.zmatLower.set(config.ematUpper, bcalc);
			this.zmatUpper = new ZMatrix(state.confSpace);
			this.zmatUpper.set(config.ematLower, bcalc);
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
			Sequence seq = confSpace.seqSpace.makeUnassignedSequence();
			for (SimpleConfSpace.Position confPos : state.confSpace.positions) {
				SeqSpace.Position seqPos = confSpace.seqSpace.getPosition(confPos.resNum);
				int rc = conf[confPos.index];
				if (seqPos != null && rc != Conf.Unassigned) {
					SimpleConfSpace.ResidueConf resConf = confPos.resConfs.get(rc);
					seq.set(seqPos, resConf.template.name);
				}
			}
			return seq;
		}

		BigDecimal calcOrderHeuristic(SimpleConfSpace.Position pos) {

			// TODO: what heuristic works best here?

			BigMath m = bigMath().set(0);
			ConfIndex index = makeConfIndex();
			for (int rc : rcs.get(pos.index)) {
				index.assignInPlace(pos.index, rc);
				m.add(calcZPathTailBound(zmatUpper, MathTools.Optimizer.Maximize, index, rcs));
				index.unassignInPlace(pos.index);
			}
			return m.div(rcs.get(pos.index).length).get();
		}

		boolean isSingleSequence(ConfIndex index) {

			for (int i=0; i<index.numUndefined; i++) {
				int pos = index.undefinedPos[i];

				// is there more than one RT at this pos?
				int rtRef = rtsByRcByPos[pos][0];
				for (int rc : rcs.get(pos)) {
					int rt = rtsByRcByPos[pos][rc];
					if (rt != rtRef) {

						// yup, multi-sequence
						return false;
					}
				}
			}

			// single-sequence
			return true;
		}

		/** WARNING: naive brute force method, for testing small trees only */
		Map<Sequence,BigInteger> countLeavesBySequence(ConfIndex index) {

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

		BigIntegerBounds boundLeavesPerSequence(ConfIndex index, RCs rcs) {

			BigIntegerBounds count = new BigIntegerBounds(BigInteger.ONE, BigInteger.ONE);

			for (int i=0; i<index.numUndefined; i++) {
				int pos = index.undefinedPos[i];

				int minCount = Integer.MAX_VALUE;
				int maxCount = 0;

				int numRts = numRtsByPos[pos];
				if (numRts > 0) {

					// count the RCs by RT
					int[] counts = new int[numRts];
					for (int rc : rcs.get(pos)) {
						int rtCount = ++counts[rtsByRcByPos[pos][rc]];
						minCount = Math.min(minCount, rtCount);
						maxCount = Math.max(maxCount, rtCount);
					}

				} else {

					// all RCs are the same RT
					minCount = rcs.getNum(pos);
					maxCount = minCount;
				}

				count.lower = count.lower.multiply(BigInteger.valueOf(minCount));
				count.upper = count.upper.multiply(BigInteger.valueOf(maxCount));
			}

			return count;
		}

		BigDecimalBounds calcZSumBounds(ConfIndex index, RCs rcs) {

			BigDecimalBounds zPathHeadBounds = calcZPathHeadBounds(index);

			if (isSingleSequence(index)) {
				return new BigDecimalBounds(
					calcZSumLower(index, rcs), // TODO: this is a really crappy lower bound! any way to make it better without minimizing?
					bigMath()
						.set(zPathHeadBounds.upper)
						.mult(calcZPathTailBound(zmatUpper, MathTools.Optimizer.Maximize, index, rcs))
						.mult(countLeafNodes(index, rcs))
						.get()
				);
			} else {
				BigIntegerBounds numLeaves = boundLeavesPerSequence(index, rcs);
				BigDecimalBounds zPathTailBounds = calcZPathTailBounds(index, rcs);
				return new BigDecimalBounds(
					bigMath() // TODO: this is a really crappy lower bound! any way to make it better?
						.set(zPathHeadBounds.lower)
						.mult(zPathTailBounds.lower)
						.mult(numLeaves.lower)
						.get(),
					bigMath()
						.set(zPathHeadBounds.upper)
						.mult(zPathTailBounds.upper)
						.mult(numLeaves.upper)
						.get()
				);
			}
		}

		BigDecimalBounds calcZSumBoundsTighter(ConfIndex index, RCs rcs, ConfDB.ConfTable confTable) {
			return calcZSumBoundsTighter(index, rcs, confTable, calcZPathHeadBoundsTighter(index, confTable));
		}

		BigDecimalBounds calcZSumBoundsTighter(ConfIndex index, RCs rcs, ConfDB.ConfTable confTable, BigDecimalBounds zPathHeadBoundsTighter) {

			// don't bother for leaf nodes
			if (index.isFullyDefined()) {
				throw new IllegalArgumentException("don't bother doing this for leaf nodes");
			}

			if (isSingleSequence(index)) {
				return new BigDecimalBounds(
					calcZSumLowerTighter(index, rcs, confTable), // minimizes
					bigMath()
						.set(zPathHeadBoundsTighter.upper)
						.mult(calcZPathTailBound(zmatUpper, MathTools.Optimizer.Maximize, index, rcs))
						.mult(countLeafNodes(index, rcs))
						.get()
				);
			} else {
				BigIntegerBounds numLeaves = boundLeavesPerSequence(index, rcs);
				BigDecimalBounds zPathTailBounds = calcZPathTailBounds(index, rcs);
				return new BigDecimalBounds(
					bigMath() // TODO: this is a really crappy lower bound! any way to make it better?
						.set(calcZPathHeadBound(zmatLower, index))
						.mult(zPathTailBounds.lower)
						.mult(numLeaves.lower)
						.get(),
					bigMath()
						.set(zPathHeadBoundsTighter.upper)
						.mult(zPathTailBounds.upper)
						.mult(numLeaves.upper)
						.get()
				);
			}
		}

		BigDecimalBounds calcZPathHeadBounds(ConfIndex index) {
			return new BigDecimalBounds(
				calcZPathHeadBound(zmatLower, index),
				calcZPathHeadBound(zmatUpper, index)
			);
		}

		BigDecimal calcZPathHeadBound(ZMatrix zmat, ConfIndex index) {

			BigMath z = bigMath().set(1.0);

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

			return z.get();
		}

		BigDecimalBounds calcZPathNodeBounds(ConfIndex index, int pos1, int rc1) {
			return new BigDecimalBounds(
				calcZPathNodeBound(zmatLower, index, pos1, rc1),
				calcZPathNodeBound(zmatUpper, index, pos1, rc1)
			);
		}

		BigDecimal calcZPathNodeBound(ZMatrix zmat, ConfIndex index, int pos1, int rc1) {

			BigMath z = bigMath().set(zmat.getOneBody(pos1, rc1));

			for (int i=0; i<index.numDefined; i++) {
				int pos2 = index.definedPos[i];
				int rc2 = index.definedRCs[i];

				z.mult(zmat.getPairwise(pos1, rc1, pos2, rc2));
			}

			return z.get();
		}

		BigDecimalBounds calcZPathTailBounds(ConfIndex index, RCs rcs) {
			return new BigDecimalBounds(
				calcZPathTailBound(zmatLower, MathTools.Optimizer.Minimize, index, rcs),
				calcZPathTailBound(zmatUpper, MathTools.Optimizer.Maximize, index, rcs)
			);
		}

		BigDecimal calcZPathTailBound(ZMatrix zmat, MathTools.Optimizer opt, ConfIndex index, RCs rcs) {

			// this is the usual A* heuristic

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

		BigDecimalBounds calcZPathBounds(ConfIndex index, RCs rcs) {
			return multBounds(
				calcZPathHeadBounds(index),
				calcZPathTailBounds(index, rcs)
			);
		}

		BigInteger countLeafNodes(ConfIndex index, RCs rcs) {

			BigInteger count = BigInteger.ONE;

			for (int i=0; i<index.numUndefined; i++) {
				int pos = index.undefinedPos[i];
				count = count.multiply(BigInteger.valueOf(rcs.getNum(pos)));
			}

			return count;
		}

		void greedilyAssignConf(ConfIndex index, RCs rcs, ZMatrix zmat) {

			// get to any leaf node and compute the subtree part of its zPath

			int numUnassigned = index.numUndefined;

			// assign each unassigned position greedily
			for (int i=0; i<numUnassigned; i++) {

				int pos = posPermutation[index.numDefined];

				// find the RC with the biggest zPathComponent
				int bestRC = -1;
				BigDecimal bestZPathNode = MathTools.BigNegativeInfinity;
				for (int rc : rcs.get(pos)) {

					BigDecimal zPathNode = calcZPathNodeBound(zmat, index, pos, rc);
					if (MathTools.isGreaterThan(zPathNode, bestZPathNode)) {
						bestZPathNode = zPathNode;
						bestRC = rc;
					}
				}

				// make the assignment
				index.assignInPlace(pos, bestRC);
			}
		}

		void unassignConf(ConfIndex index, int numAssigned) {

			// undo all the assignments in the reverse order they were assigned
			for (int i=index.numPos-1; i>=numAssigned; i--) {
				index.unassignInPlace(posPermutation[i]);
			}
		}

		BigDecimal calcZSumLower(ConfIndex index, RCs rcs) {

			// compute the zPathLower for an arbitrary leaf node
			// but heuristically try to pick one with a high zPathLower
			int numAssigned = index.numDefined;
			greedilyAssignConf(index, rcs, zmatLower);
			BigDecimal zPathLower = calcZPathHeadBound(zmatLower, index);
			unassignConf(index, numAssigned);
			return zPathLower;
		}

		BigDecimalBounds calcZPathHeadBoundsTighter(ConfIndex index, ConfDB.ConfTable confTable) {
			return new BigDecimalBounds(
				calcZPathHeadBound(zmatLower, index),
				calcZPathHeadUpperTighter(index, confTable)
			);
		}

		BigDecimal calcZPathHeadUpperTighter(ConfIndex index, ConfDB.ConfTable confTable) {
			RCTuple tuple = new RCTuple(index);
			ResidueInteractions inters = confEcalc.makeTupleInters(tuple);
			double e = confEcalc.calcEnergy(tuple, inters, confTable);
			return bcalc.calcPrecise(e);
		}

		BigDecimal calcZPath(ConfIndex index, ConfDB.ConfTable confTable) {

			if (!index.isFullyDefined()) {
				throw new IllegalArgumentException("not a full conf");
			}

			double e = confEcalc.calcEnergy(new RCTuple(index), confTable);
			// TEMP
			//return bcalc.calcPrecise(e);
			BigDecimal z = bcalc.calcPrecise(e);
			log("%s  conf %4d  E = %12.6f, Z = %.6e", index, minimizationCounter.incrementAndGet(), e, z);

			return z;
		}

		/** a reminder to myself this this is a bad idea */
		@Deprecated
		BigDecimal calcZPathHead(ConfIndex index, ConfDB.ConfTable confTable) {
			throw new Error("stop trying to do this! There's no such thing!");
		}


		BigDecimal calcZSumLowerTighter(ConfIndex index, RCs rcs, ConfDB.ConfTable confTable) {

			// compute the zPath for an arbitrary leaf node
			// but heuristically try to pick one with a high zPath
			int numAssigned = index.numDefined;
			greedilyAssignConf(index, rcs, zmatUpper);
			BigDecimal zPath = calcZPath(index, confTable);
			unassignConf(index, numAssigned);
			return zPath;
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

	// TEMP
	private AtomicInteger minimizationCounter = new AtomicInteger(0);

	private BigDecimalBounds multBounds(BigDecimalBounds a, BigDecimalBounds b) {
		return new BigDecimalBounds(
			bigMath().set(a.lower).mult(b.lower).get(),
			bigMath().set(a.upper).mult(b.upper).get()
		);
	}
}
