package edu.duke.cs.osprey.sofea;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.lute.LUTEConfEnergyCalculator;
import edu.duke.cs.osprey.pruning.PruningMatrix;
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
 * TODO: doc me!
 */
public class Sofea {

	public static class Builder {

		public final MultiStateConfSpace confSpace;
		public final Criterion criterion;

		private StateConfig[] stateConfigs;
		private MathContext mathContext = new MathContext(32, RoundingMode.HALF_UP);
		private File seqdbFile = new File("seq.db");
		private MathContext seqdbMathContext = new MathContext(128, RoundingMode.HALF_UP);
		private File fringedbFile = new File("fringe.db");
		private int fringedbMiB = 10;
		private boolean showProgress = true;

		// NOTE: don't need much precision for most math, but need lots of precision for seqdb math

		public Builder(MultiStateConfSpace confSpace, Criterion criterion) {

			this.confSpace = confSpace;
			this.criterion = criterion;

			this.stateConfigs = new StateConfig[confSpace.states.size()];
		}

		public Builder configEachState(Function<MultiStateConfSpace.State,StateConfig> configurator) {
			for (MultiStateConfSpace.State state : confSpace.states) {
				stateConfigs[state.index] = configurator.apply(state);
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

		public Builder setFringeDBMiB(int val) {
			fringedbMiB = val;
			return this;
		}

		public Builder setShowProgress(boolean val) {
			showProgress = val;
			return this;
		}

		public Sofea make() {

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
				criterion,
				Arrays.asList(stateConfigs),
				mathContext,
				seqdbFile,
				seqdbMathContext,
				fringedbFile,
				fringedbMiB,
				showProgress
			);
		}
	}

	public static class StateConfig {

		public final LUTEConfEnergyCalculator luteEcalc;
		public final PruningMatrix pmat;

		public StateConfig(LUTEConfEnergyCalculator luteEcalc, PruningMatrix pmat) {
			this.luteEcalc = luteEcalc;
			this.pmat = pmat;
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

	/** decides if computation should continue or not */
	public static interface Criterion {
		boolean isFinished(SeqDB seqdb, FringeDB fringedb, long sweepCount);
	}


	public final MultiStateConfSpace confSpace;
	public final Criterion criterion;
	public final List<StateConfig> stateConfigs;
	public final MathContext mathContext;
	public final File seqdbFile;
	public final MathContext seqdbMathContext;
	public final File fringedbFile;
	public final int fringedbMiB;
	public final boolean showProgress;

	private final List<StateInfo> stateInfos;

	private Sofea(MultiStateConfSpace confSpace, Criterion criterion, List<StateConfig> stateConfigs, MathContext mathContext, File seqdbFile, MathContext seqdbMathContext, File fringedbFile, int fringedbMiB, boolean showProgress) {

		this.confSpace = confSpace;
		this.criterion = criterion;
		this.stateConfigs = stateConfigs;
		this.mathContext = mathContext;
		this.seqdbFile = seqdbFile;
		this.seqdbMathContext = seqdbMathContext;
		this.fringedbFile = fringedbFile;
		this.fringedbMiB = fringedbMiB;
		this.showProgress = showProgress;

		// TODO: support single tuples for conf spaces
		// TEMP: blow up if we get single tuples
		for (MultiStateConfSpace.State state : confSpace.states) {
			if (state.confSpace.positions.size() <= 1) {
				throw new Error("TODO: support small conf spaces");
			}
		}

		// init the state info
		stateInfos = confSpace.states.stream()
			.map(state -> new StateInfo(state))
			.collect(Collectors.toList());
	}

	public StateConfig getConfig(MultiStateConfSpace.State state) {
		return stateConfigs.get(state.index);
	}

	private BigMath bigMath() {
		return new BigMath(mathContext);
	}

	public SeqDB openSeqDB() {
		return new SeqDB(confSpace, seqdbMathContext, seqdbFile);
	}

	public FringeDB openFringeDB() {
		if (fringedbFile.exists()) {
			return FringeDB.open(confSpace, fringedbFile);
		} else {
			return FringeDB.create(confSpace, fringedbFile, fringedbMiB*1024*1024, mathContext);
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

		// OPTIMIZATION: allocate one RCTuple instance for all the triple lookups
		// and set positions in 3-2-1 order so the positions are sorted in ascending order
		RCTuple tripleTuple = new RCTuple(0, 0, 0, 0, 0, 0);

		try (SeqDB seqdb = openSeqDB()) {
			try (FringeDB fringedb = openFringeDB()) {

				// process the root node for each state
				FringeDB.Transaction fringetx = fringedb.transaction();
				SeqDB.Transaction seqtx = seqdb.transaction();
				for (MultiStateConfSpace.State state : confSpace.states) {
					StateInfo stateInfo = stateInfos.get(state.index);

					// get a multi-sequence Z bound on the root node
					MathTools.BigDecimalBounds rootBound;
					{
						ConfIndex index = stateInfo.makeConfIndex();
						BigDecimal minZ = stateInfo.optimizeZ(index, tripleTuple, MathTools.Optimizer.Minimize);
						BigDecimal maxZ = stateInfo.optimizeZ(index, tripleTuple, MathTools.Optimizer.Maximize);
						MathTools.BigIntegerBounds size = stateInfo.count(index);

						rootBound = new MathTools.BigDecimalBounds(
							bigMath()
								.set(minZ)
								.mult(size.lower)
								.mult(stateInfo.blute.factor)
								.get(),
							bigMath()
								.set(maxZ)
								.mult(size.upper)
								.mult(stateInfo.blute.factor)
								.get()
						);
					}

					// init the fringe with the root node
					fringetx.writeRootNode(state, rootBound, stateInfo.blute.factor);
					seqtx.addZ(state, state.confSpace.makeUnassignedSequence(), rootBound);
				}
				fringetx.commit();
				fringedb.finishSweep();
				seqtx.commit();
			}
		}
	}

	private static class Node {

		final int[] conf;
		final BigDecimalBounds zbounds;
		final BigDecimal zpath;

		Node(int[] conf, BigDecimalBounds zbounds, BigDecimal zpath) {
			this.conf = conf;
			this.zbounds = zbounds;
			this.zpath = zpath;
		}
	}

	private static class ZVal {

		final int[] conf;
		final BigDecimal z;

		ZVal(int[] conf, BigDecimal z) {
			this.conf = conf;
			this.z = z;
		}
	}

	private class NodeTransaction {

		final FringeDB.Transaction fringetx;
		final SeqDB.Transaction seqtx;
		final MultiStateConfSpace.State state;

		final List<Node> replacementNodes = new ArrayList<>();
		final List<ZVal> zvals = new ArrayList<>();

		NodeTransaction(FringeDB.Transaction fringetx, SeqDB.Transaction seqtx, MultiStateConfSpace.State state) {
			this.fringetx = fringetx;
			this.seqtx = seqtx;
			this.state = state;
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

		void addReplacementNode(Node node) {
			replacementNodes.add(node);
		}

		void addReplacementNode(ConfIndex index, BigDecimalBounds zbounds, BigDecimal zpath) {
			replacementNodes.add(new Node(Conf.make(index), zbounds, zpath));
		}

		int numReplacementNodes() {
			return replacementNodes.size();
		}

		void addZVal(ConfIndex index, BigDecimal z) {
			zvals.add(new ZVal(Conf.make(index), z));
		}

		boolean hasRoomToReplace() {
			return fringetx.dbHasRoomFor(replacementNodes.size());
		}

		void replace(Node node) {

			flushTransactionsIfNeeded();

			// subtract the node zbounds
			seqtx.subZ(state, makeSeq(node.conf), node.zbounds);

			// update fringedb and seqdb with the replacement nodes
			for (Node replacementNode : replacementNodes) {
				fringetx.writeReplacementNode(replacementNode.conf, replacementNode.zbounds, replacementNode.zpath);
				seqtx.addZ(state, makeSeq(replacementNode.conf), replacementNode.zbounds);
			}

			// add the z values
			for (ZVal zval : zvals) {
				seqtx.addZ(state, makeSeq(zval.conf), zval.z);
			}
		}

		void requeue(Node node) {

			flushTransactionsIfNeeded();

			// ignore all of the seqdb changes

			// move the node to the end of the fringedb queue
			fringetx.writeReplacementNode(node.conf, node.zbounds, node.zpath);
		}

		private void flushTransactionsIfNeeded() {

			// skip if the fringe db transaction isn't full
			if (fringetx.txHasRoomFor(replacementNodes.size())) {
				return;
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
		}
	}

	private static class StateStats {
		long read = 0;
		long expanded = 0;
		long replaced = 0;
		long added = 0;
		long requeued = 0;
	}

	public void refine() {
		try (SeqDB seqdb = openSeqDB()) {
			try (FringeDB fringedb = openFringeDB()) {

				// get the initial zmax
				BigDecimal[] zmax = confSpace.states.stream()
					.map(state -> fringedb.getZMax(state))
					.toArray(size -> new BigDecimal[size]);

				final double zmaxFactor = Math.pow(Math.E, 4.0);

				// OPTIMIZATION: allocate one RCTuple instance for all the triple lookups
				// and set positions in 3-2-1 order so the positions are sorted in ascending order
				RCTuple tripleTuple = new RCTuple(0, 0, 0, 0, 0, 0);

				long sweepCount = 0;
				while (true) {

					// check the termination criterion
					if (criterion != null && criterion.isFinished(seqdb, fringedb, sweepCount)) {
						log("SOFEA finished, criterion satisfied");
						break;
					}

					// stop if we ran out of fringe nodes
					if (fringedb.isEmpty()) {
						log("SOFEA finished, explored every node");
						break;
					}

					// reduce zmax before each sweep
					for (MultiStateConfSpace.State state : confSpace.states) {
						zmax[state.index] = bigMath()
							.set(zmax[state.index])
							.div(zmaxFactor)
							.get();
					}

					// show progress if needed
					sweepCount++;
					if (showProgress) {
						log("sweep %d", sweepCount);
					}

					// init sweep stats
					StateStats[] stats = new StateStats[confSpace.states.size()];
					for (MultiStateConfSpace.State state : confSpace.states) {
						stats[state.index] = new StateStats();
					}

					// start (or resume) the sweep
					FringeDB.Transaction fringetx = fringedb.transaction();
					SeqDB.Transaction seqtx = seqdb.transaction();
					Progress progress = new Progress(fringetx.numNodesToRead());
					while (fringetx.hasNodesToRead()) {

						// read the next node
						fringetx.readNode();
						MultiStateConfSpace.State state = fringetx.state();
						Node node = new Node(
							fringetx.conf(),
							fringetx.zbounds(),
							fringetx.zpath()
						);

						NodeTransaction nodetx = new NodeTransaction(fringetx, seqtx, state);
						stats[nodetx.state.index].read++;

						ConfIndex index = new ConfIndex(state.confSpace.positions.size());

						boolean wasExpanded = design(nodetx, zmax[nodetx.state.index], node, tripleTuple, index);
						if (wasExpanded) {
							stats[nodetx.state.index].expanded++;
						}

						if (nodetx.hasRoomToReplace()) {
							stats[nodetx.state.index].added += nodetx.numReplacementNodes();
							stats[nodetx.state.index].replaced++;
							nodetx.replace(node);
						} else {
							stats[nodetx.state.index].requeued++;
							nodetx.requeue(node);
						}

						if (showProgress) {
							progress.incrementProgress();
						}
					}
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
			}
		}
	}

	private boolean design(NodeTransaction nodetx, BigDecimal zmax, Node node, RCTuple tripleTuple, ConfIndex index) {

		// skip this tree if it's too small, and add the node to the fringe set
		if (MathTools.isLessThan(node.zbounds.upper, zmax)) {
			nodetx.addReplacementNode(node);
			return false;
		}

		// otherwise, recurse

		// update the conf index
		Conf.index(node.conf, index);

		StateInfo stateInfo = stateInfos.get(nodetx.state.index);
		boolean isRoot = index.numDefined == 0;
		boolean isRCLeaf = index.numDefined + 1 == index.numPos;

		int pos = index.numDefined;
		for (int rc : stateInfo.rcs.get(pos)) {

			// update the zpath with this RC
			BigDecimal zpathrc = node.zpath;
			if (!isRoot) {

				BigDecimal zrc = stateInfo.getZPart(index, tripleTuple, pos, rc);

				// this subtree contributes nothing to Z
				if (MathTools.isZero(zrc)) {
					continue;
				}

				zpathrc = bigMath()
					.set(zpathrc)
					.mult(zrc)
					.get();
			}

			if (isRCLeaf) {

				assert (!isRoot);

				// hit bottom, add the leaf node to the seqdb
				index.assignInPlace(pos, rc);
				nodetx.addZVal(index, zpathrc);
				index.unassignInPlace(pos);

			} else {

				index.assignInPlace(pos, rc);

				// get the subtree bounds
				BigDecimal minZ = stateInfo.optimizeZ(index, tripleTuple, MathTools.Optimizer.Minimize);
				BigDecimal maxZ = stateInfo.optimizeZ(index, tripleTuple, MathTools.Optimizer.Maximize);
				MathTools.BigIntegerBounds count = stateInfo.count(index);

				// if the upper bound is zero, this subtree contributes nothing to Z
				if (MathTools.isZero(maxZ) || count.upper.compareTo(BigInteger.ZERO) == 0) {
					index.unassignInPlace(pos);
					continue;
				}

				MathTools.BigDecimalBounds boundsrc  = new MathTools.BigDecimalBounds(
					bigMath()
						.set(minZ)
						.mult(count.lower)
						.mult(zpathrc)
						.get(),
					bigMath()
						.set(maxZ)
						.mult(count.upper)
						.mult(zpathrc)
						.get()
				);

				// recurse
				design(nodetx, zmax, new Node(Conf.make(index), boundsrc, zpathrc), tripleTuple, index);

				index.unassignInPlace(pos);
			}
		}

		return true;
	}

	public BigDecimal calcZ(MultiStateConfSpace.State state, Sequence seq) {
		return stateInfos.get(state.index).calcZ(seq);
	}

	private class StateInfo {

		final MultiStateConfSpace.State state;
		final BoltzmannLute blute;
		final PruningMatrix pmat;
		final RCs rcs;
		final List<SimpleConfSpace.Position> positions;

		private final TupleMatrixGeneric<BigDecimal> optrc3Energies;
		private final int[][] rtsByRcByPos;
		private final int[] numRtsByPos;

		StateInfo(MultiStateConfSpace.State state) {

			this.state = state;
			StateConfig config = stateConfigs.get(state.index);
			this.blute = new BoltzmannLute(config.luteEcalc, mathContext);
			this.pmat = config.pmat;
			this.rcs = new RCs(state.confSpace);
			this.positions = state.confSpace.positions;

			// OPTIMIZATION: allocate one RCTuple instance for all the triple lookups
			// and set positions in 3-2-1 order so the positions are sorted in ascending order
			RCTuple tripleTuple = new RCTuple(0, 0, 0, 0, 0, 0);

			// pre-calculate all the optimal triple Z values for each RC pair
			optrc3Energies = new TupleMatrixGeneric<>(state.confSpace);
			for (int pos1=0; pos1<rcs.getNumPos(); pos1++) {
				for (int rc1 : rcs.get(pos1)) {
					for (int pos2=0; pos2<pos1; pos2++) {
						for (int rc2 : rcs.get(pos2)) {

							BigDecimal optrc3Energy = MathTools.Optimizer.Maximize.initBigDecimal();
							for (int pos3=0; pos3<pos2; pos3++) {
								for (int rc3 : rcs.get(pos3)) {

									tripleTuple.set(pos3, rc3, pos2, rc2, pos1, rc1);
									BigDecimal triple = blute.get(tripleTuple);
									if (triple != null) {
										optrc3Energy = MathTools.Optimizer.Maximize.opt(optrc3Energy, triple);
									}
								}
							}

							optrc3Energies.setPairwise(pos1, rc1, pos2, rc2, optrc3Energy);
						}
					}
				}
			}

			// calculate all the RTs by RC and pos
			rtsByRcByPos = new int[positions.size()][];
			numRtsByPos = new int[positions.size()];
			for (SimpleConfSpace.Position confPos : positions) {
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
			ConfIndex index = new ConfIndex(positions.size());
			index.updateUndefined();
			return index;
		}

		BigDecimal getZPart(ConfIndex confIndex, RCTuple tripleTuple, int pos1, int rc1) {

			assert (confIndex.numDefined > 0);

			if (pmat.getOneBody(pos1, rc1)) {
				return BigDecimal.ZERO;
			}

			BigMath math = bigMath().set(1.0);
			for (int i=0; i<confIndex.numDefined; i++) {
				int pos2 = confIndex.definedPos[i];
				int rc2 = confIndex.definedRCs[i];

				if (pmat.getPairwise(pos1, rc1, pos2, rc2)) {
					return BigDecimal.ZERO;
				}

				math.mult(blute.get(pos1, rc1, pos2, rc2));

				for (int j=0; j<i; j++) {
					int pos3 = confIndex.definedPos[j];
					int rc3 = confIndex.definedRCs[j];

					if (pmat.getPairwise(pos1, rc1, pos3, rc3)) {
						return BigDecimal.ZERO;
					}
					if (pmat.getPairwise(pos2, rc2, pos3, rc3)) {
						return BigDecimal.ZERO;
					}
					tripleTuple.set(pos3, rc3, pos2, rc2, pos1, rc1);
					if (pmat.getTuple(tripleTuple)) {
						return BigDecimal.ZERO;
					}

					BigDecimal triple = blute.get(tripleTuple);
					if (triple != null) {
						math.mult(triple);
					}
				}
			}

			return math.get();
		}

		BigDecimal optimizeZ(ConfIndex index, RCTuple tripleTuple, MathTools.Optimizer opt) {

			BigMath energy = bigMath().set(1.0);

			// get the score for each undefined position
			for (int i=0; i<index.numUndefined; i++) {
				int pos1 = index.undefinedPos[i];

				// optimize over possible assignments to pos1
				BigDecimal pos1Energy = opt.initBigDecimal();
				for (int rc1 : rcs.get(pos1)) {

					BigMath rc1Energy = bigMath();
					if (isPruned(index, tripleTuple, pos1, rc1)) {

						// prune tuple, no contribution to Z
						rc1Energy.set(0.0);

					} else {

						rc1Energy.set(1.0);

						// interactions with defined residues
						for (int j=0; j<index.numDefined; j++) {
							int pos2 = index.definedPos[j];
							int rc2 = index.definedRCs[j];

							rc1Energy.mult(blute.get(pos1, rc1, pos2, rc2));

							for (int k=0; k<j; k++) {
								int pos3 = index.definedPos[k];
								int rc3 = index.definedRCs[k];

								// triples are optional
								tripleTuple.set(pos3, rc3, pos2, rc2, pos1, rc1);
								BigDecimal triple = blute.get(tripleTuple);
								if (triple != null) {
									rc1Energy.mult(triple);
								}
							}
						}

						// interactions with undefined residues
						for (int j=0; j<i; j++) {
							int pos2 = index.undefinedPos[j];

							// optimize over possible assignments to pos2
							BigDecimal optrc2Energy = opt.initBigDecimal();
							for (int rc2 : rcs.get(pos2)) {

								BigMath rc2Energy = bigMath();

								if (pmat.getPairwise(pos1, rc1, pos2, rc2)) {

									// prune tuple, no contribution to Z
									rc2Energy.set(0.0);

								} else {

									// pair with pos2
									rc2Energy.set(blute.get(pos1, rc1, pos2, rc2));

									// triples with defined positions
									for (int k=0; k<index.numDefined; k++) {
										int pos3 = index.definedPos[k];
										int rc3 = index.definedRCs[k];

										// triples are optional
										tripleTuple.set(pos3, rc3, pos2, rc2, pos1, rc1);
										BigDecimal triple = blute.get(tripleTuple);
										if (triple != null) {
											rc2Energy.mult(triple);
										}
									}

									// triples with undefined positions
									for (int k=0; k<j; k++) {
										int pos3 = index.undefinedPos[k];
										BigDecimal optrc3Energy = optrc3Energies.getPairwise(pos1, rc1, pos2, rc2);
										if (MathTools.isFinite(optrc3Energy)) {
											rc2Energy.mult(optrc3Energy);
										}
									}
								}

								optrc2Energy = opt.opt(optrc2Energy, rc2Energy.get());
							}

							rc1Energy.mult(optrc2Energy);
						}
					}

					pos1Energy = opt.opt(pos1Energy, rc1Energy.get());
				}

				assert (MathTools.isFinite(pos1Energy));
				energy.mult(pos1Energy);
			}

			return energy.get();
		}

		boolean isPruned(ConfIndex index, RCTuple tripleTuple, int pos1, int rc1) {

			if (pmat.getOneBody(pos1, rc1)) {
				return true;
			}

			for (int i=0; i<index.numDefined; i++) {
				int pos2 = index.definedPos[i];
				int rc2 = index.definedRCs[i];

				if (pmat.getPairwise(pos1, rc1, pos2, rc2)) {
					return true;
				}

				for (int j=0; j<i; j++) {
					int pos3 = index.definedPos[j];
					int rc3 = index.definedRCs[j];

					if (pmat.getPairwise(pos1, rc1, pos2, rc2)) {
						return true;
					}
					if (pmat.getPairwise(pos1, rc1, pos3, rc3)) {
						return true;
					}
					tripleTuple.set(pos3, rc3, pos2, rc2, pos1, rc1);
					if (pmat.getTuple(tripleTuple)) {
						return true;
					}
				}
			}

			return false;
		}

		MathTools.BigIntegerBounds count(ConfIndex index) {

			MathTools.BigIntegerBounds count = new MathTools.BigIntegerBounds(BigInteger.ONE, BigInteger.ONE);

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

		BigDecimal calcZ() {
			ConfIndex index = makeConfIndex();
			return calcZ(index, rcs);
		}

		BigDecimal calcZ(Sequence seq) {
			ConfIndex index = makeConfIndex();
			RCs rcs = seq.makeRCs(state.confSpace);
			return calcZ(index, rcs);
		}

		BigDecimal calcZ(ConfIndex index, RCs rcs) {

			// cannot work on leaf nodes
			assert (index.numDefined < index.numPos);

			BigDecimal z = BigDecimal.ZERO;

			// OPTIMIZATION: allocate one RCTuple instance for all the triple lookups
			// and set positions in 3-2-1 order so the positions are sorted in ascending order
			RCTuple tripleTuple = new RCTuple(0, 0, 0, 0, 0, 0);

			// pick the next pos to assign
			int pos = index.numDefined;

			for (int rc : rcs.get(pos)) {

				// get the zpart
				BigDecimal zpart;
				if (index.numDefined == 0) {
					zpart = blute.factor;
				} else {
					zpart = getZPart(index, tripleTuple, pos, rc);
				}

				// short circuit for efficiency
				if (MathTools.isZero(zpart)) {
					continue;
				}

				if (index.numDefined < index.numPos - 1) {

					// have positions left to assign, so recurse
					index.assignInPlace(pos, rc);
					z = bigMath()
						.set(calcZ(index, rcs))
						.mult(zpart)
						.add(z)
						.get();
					index.unassignInPlace(pos);

				} else {

					// all positions assigned, just add the zpart
					z = bigMath()
						.set(zpart)
						.add(z)
						.get();
				}
			}

			return z;
		}
	}
}
