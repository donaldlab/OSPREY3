package edu.duke.cs.osprey.sofea;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.lute.LUTEConfEnergyCalculator;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.Log;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;
import edu.duke.cs.osprey.tools.resultdoc.Plot;
import edu.duke.cs.osprey.tools.resultdoc.ResultDoc;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;
import static edu.duke.cs.osprey.tools.Log.logf;


/**
 * TODO: doc me!
 */
public class Sofea {

	public static class Builder {

		public final MultiStateConfSpace confSpace;
		public final MultiStateConfSpace.LMFE objective;

		private StateConfig[] stateConfigs;
		private MathContext mathContext = new MathContext(32, RoundingMode.HALF_UP);
		private File seqdbFile = new File("seq.db");
		private MathContext seqdbMathContext = new MathContext(256, RoundingMode.HALF_UP);
		private File fringedbFile = new File("fringe.db");
		private int fringedbMiB = 10;

		// NOTE: don't need much precision for most math, but need lots of precision for seqdb math

		public Builder(MultiStateConfSpace confSpace, MultiStateConfSpace.LMFE objective) {

			this.confSpace = confSpace;
			this.objective = objective;

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
				objective,
				Arrays.asList(stateConfigs),
				mathContext,
				seqdbFile,
				seqdbMathContext,
				fringedbFile,
				fringedbMiB
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


	public final MultiStateConfSpace confSpace;
	public final MultiStateConfSpace.LMFE objective;
	public final List<StateConfig> stateConfigs;
	public final MathContext mathContext;
	public final File seqdbFile;
	public final MathContext seqdbMathContext;
	public final File fringedbFile;
	public final int fringedbMiB;

	private final List<StateInfo> stateInfos;

	private Sofea(MultiStateConfSpace confSpace, MultiStateConfSpace.LMFE objective, List<StateConfig> stateConfigs, MathContext mathContext, File seqdbFile, MathContext seqdbMathContext, File fringedbFile, int fringedbMiB) {

		this.confSpace = confSpace;
		this.objective = objective;
		this.stateConfigs = stateConfigs;
		this.mathContext = mathContext;
		this.seqdbFile = seqdbFile;
		this.seqdbMathContext = seqdbMathContext;
		this.fringedbFile = fringedbFile;
		this.fringedbMiB = fringedbMiB;

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

		try (SeqDB seqdb = openSeqDB()) {
			try (FringeDB fringedb = openFringeDB()) {

				// process the root node for each state
				FringeDB.Roots roots = fringedb.roots();
				SeqDB.Transaction seqtx = seqdb.transaction();
				for (MultiStateConfSpace.State state : confSpace.states) {
					StateInfo stateInfo = stateInfos.get(state.index);

					// get a multi-sequence Z bound on the root node
					MathTools.BigDecimalBounds rootBound;
					{
						ConfIndex index = stateInfo.makeConfIndex();
						BigDecimal minZ = stateInfo.optimizeZ(index, MathTools.Optimizer.Minimize);
						BigDecimal maxZ = stateInfo.optimizeZ(index, MathTools.Optimizer.Maximize);
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
					roots.set(state, rootBound, stateInfo.blute.factor);
					ConfIndex index = stateInfo.makeConfIndex();
					seqtx.addZ(state, index, rootBound);
				}
				roots.commit();
				seqtx.commit();

				// TEMP
				log("seed");
				log("\tfringe size: %d/%d", fringedb.getNumNodes(), fringedb.getCapacity());
			}
		}
	}

	public void refine() {

		// TODO: add stopping criteria of some kind?

		try (SeqDB seqdb = openSeqDB()) {
			try (FringeDB fringedb = openFringeDB()) {

				final double zmaxFactor = Math.pow(Math.E, 4.0);

				// TEMP: try a few operations
				for (int i=0; i<8; i++) {

					// reduce zmax each iteration
					List<BigDecimal> zmax = confSpace.states.stream()
						.map(state -> {
							BigDecimal stateZMax = fringedb.getZMax(state);
							if (stateZMax == null) {
								return null;
							} else {
								return bigMath()
									.set(fringedb.getZMax(state))
									.div(zmaxFactor)
									.get();
							}
						})
						.collect(Collectors.toList());

					// TEMP
					log("op %3d", i);
					for (MultiStateConfSpace.State state : confSpace.states) {
						log("\tzmax=%s  /f=%s", Log.formatBigLn(fringedb.getZMax(state)), Log.formatBigLn(zmax.get(state.index)));
					}

					long[] numExpanded = { 0 };
					long[] numRead = { 0 };

					FringeDB.Sweep sweep = fringedb.sweep();
					SeqDB.Transaction seqtx = seqdb.transaction();
					while (!sweep.isEmpty()) {

						sweep.read();
						MultiStateConfSpace.State state = sweep.state();
						ConfIndex index = sweep.index();
						BigDecimalBounds bounds = sweep.bounds();
						BigDecimal zpath = sweep.zpath();

						seqtx.subZ(state, index, bounds);

						boolean wasExpanded = design(state, index, bounds, zmax.get(state.index), zpath, sweep, seqtx);
						if (wasExpanded) {
							numExpanded[0]++;
						}
						numRead[0]++;

						if (sweep.hasSpaceToReplace()) {
							sweep.replace();
						} else {
							sweep.requeueAndDiscardChildren();
						}

						// TODO: come up with better batching heuristic
						if (sweep.isEmpty() || numRead[0] % 100 == 0) {
							sweep.commit();
							seqtx.commit();
						}
					}
					sweep.finish();

					// TEMP
					log("\texpanded %d/%d, fringe size: %d/%d",
						numExpanded[0], numRead[0],
						fringedb.getNumNodes(), fringedb.getCapacity()
					);

					// TEMP
					evaluateObjective(seqdb);

					// stop if we ran out of fringe nodes
					if (fringedb.isEmpty()) {
						break;
					}
				}
			}
		}
	}

	// TEMP: try to evaluate the objective function
	private void evaluateObjective(SeqDB seqdb) {

		log("objective function:");

		BoltzmannCalculator bcalc = new BoltzmannCalculator(mathContext);

		// get the unsequenced z values
		DoubleBounds[] unsequencedFreeEnergy = new DoubleBounds[confSpace.unsequencedStates.size()];
		for (MultiStateConfSpace.State state : confSpace.unsequencedStates) {
			unsequencedFreeEnergy[state.unsequencedIndex] = bcalc.freeEnergyPrecise(seqdb.getUnsequenced(state.unsequencedIndex));
		}

		// for each sequence and partial sequence encountered so far...
		for (Map.Entry<Sequence,SeqDB.SeqInfo> entry : seqdb.getSequences()) {

			Sequence seq = entry.getKey();
			SeqDB.SeqInfo seqInfo = entry.getValue();

			// compute bounds on the objective function
			DoubleBounds objectiveBounds = new DoubleBounds(0.0, 0.0);
			for (MultiStateConfSpace.State state : objective.states()) {

				DoubleBounds g;
				if (state.isSequenced) {
					g = bcalc.freeEnergyPrecise(seqInfo.bounds[state.sequencedIndex]);
				} else {
					g = unsequencedFreeEnergy[state.unsequencedIndex];
				}

				objectiveBounds.lower += g.lower*objective.getWeight(state);
				objectiveBounds.upper += g.upper*objective.getWeight(state);
			}

			logf("\tobj=%s", objectiveBounds.toString(4, 9));
			for (MultiStateConfSpace.State state : seqdb.confSpace.sequencedStates) {
				BigDecimalBounds stateBounds = seqInfo.bounds[state.sequencedIndex];

				logf("    %s=%s", state.name, Log.formatBigLn(seqInfo.bounds[state.sequencedIndex]));
				if (seq.isFullyAssigned()) {
					logf(" %.4f", stateBounds.delta(mathContext));
				}
			}
			log("    [%s]", seq);
		}
	}

	// TODO: move to reporter class of some kind?
	public void makeResultDoc(File file) {

		try (SeqDB seqdb = openSeqDB()) {

			BoltzmannCalculator bcalc = new BoltzmannCalculator(seqdbMathContext);
			try (ResultDoc doc = new ResultDoc(file)) {

				doc.h1("SOFEA Results");

				for (MultiStateConfSpace.State state : confSpace.states) {

					doc.println();
					doc.h2(state.name);

					Plot plot = new Plot();
					plot.ylabels = new ArrayList<>();

					Plot.Intervals intervals = plot.new IntervalsX();
					intervals.name = "Bounds";
					intervals.data = new ArrayList<>();

					Plot.PointsX points = plot.new PointsX();
					points.name = "Exact Value";
					points.data = new ArrayList<>();

					if (state.isSequenced) {

						for (Map.Entry<Sequence,SeqDB.SeqInfo> entry : seqdb.getSequences()) {

							Sequence seq = entry.getKey();
							plot.ylabels.add(seq.toString());

							SeqDB.SeqInfo seqInfo = entry.getValue();
							intervals.data.add(new Plot.Interval(
								bcalc.freeEnergyPrecise(seqInfo.bounds[state.sequencedIndex].upper),
								bcalc.freeEnergyPrecise(seqInfo.bounds[state.sequencedIndex].lower)
							));

							// TODO: freeEnergyPrecise is really slow for some reason

							// TEMP: show exact values
							if (seq.isFullyAssigned()) {
								BigDecimal z = stateInfos.get(state.index).calcZ(seq);
								points.data.add(bcalc.freeEnergyPrecise(z));
							} else {
								points.data.add(null);
							}
						}

					} else {

						plot.ylabels.add("");

						MathTools.BigDecimalBounds bounds = seqdb.getUnsequenced(state.unsequencedIndex);
						intervals.data.add(new Plot.Interval(
							bcalc.freeEnergyPrecise(bounds.upper),
							bcalc.freeEnergyPrecise(bounds.lower)
						));

						// TEMP: show exact values
						BigDecimal z = stateInfos.get(state.index).calcZ();
						points.data.add(bcalc.freeEnergyPrecise(z));
					}

					plot.key = "on tmargin horizontal";
					plot.xlabel = "Free Energy (kcal/mol)";
					plot.xlabelrotate = 30.0;
					plot.height = 70 + plot.ylabels.size()*40;
					doc.plot(plot);
				}
			}
		}
	}

	private boolean design(MultiStateConfSpace.State state, ConfIndex index, MathTools.BigDecimalBounds bounds, BigDecimal zmax, BigDecimal zpath, FringeDB.Sweep sweep, SeqDB.Transaction seqtx) {

		// skip this tree if it's too small, and add the node to the fringe set
		if (MathTools.isLessThan(bounds.upper, zmax)) {
			sweep.addChild(state, index, bounds, zpath);
			seqtx.addZ(state, index, bounds);
			return false;
		}

		// otherwise, recurse

		StateInfo stateInfo = stateInfos.get(state.index);

		boolean isRoot = index.numDefined == 0;
		boolean isRCLeaf = index.numDefined + 1 == index.numPos;

		int pos = index.numDefined;
		for (int rc : stateInfo.rcs.get(pos)) {

			// update the zpath with this RC
			BigDecimal zpathrc = zpath;
			if (!isRoot) {

				BigDecimal zrc = stateInfo.getZPart(index, pos, rc);

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
				// TODO: can optimize seq creation?
				index.assignInPlace(pos, rc);
				seqtx.addZ(state, index, zpathrc);
				index.unassignInPlace(pos);

			} else {

				// get the subtree bounds
				index.assignInPlace(pos, rc);
				BigDecimal minZ = stateInfo.optimizeZ(index, MathTools.Optimizer.Minimize);
				BigDecimal maxZ = stateInfo.optimizeZ(index, MathTools.Optimizer.Maximize);
				MathTools.BigIntegerBounds count = stateInfo.count(index);
				index.unassignInPlace(pos);

				// if the upper bound is zero, this subtree contributes nothing to Z
				if (MathTools.isZero(maxZ) || count.upper.compareTo(BigInteger.ZERO) == 0) {
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
				index.assignInPlace(pos, rc);
				design(state, index, boundsrc, zmax, zpathrc, sweep, seqtx);
				index.unassignInPlace(pos);
			}
		}

		return true;
	}

	private class StateInfo {

		final MultiStateConfSpace.State state;
		final BoltzmannLute blute;
		final PruningMatrix pmat;
		final RCs rcs;
		final List<SimpleConfSpace.Position> positions;

		StateInfo(MultiStateConfSpace.State state) {
			this.state = state;
			StateConfig config = stateConfigs.get(state.index);
			this.blute = new BoltzmannLute(config.luteEcalc, mathContext);
			this.pmat = config.pmat;
			this.rcs = new RCs(state.confSpace);
			this.positions = state.confSpace.positions;
		}

		ConfIndex makeConfIndex() {
			ConfIndex index = new ConfIndex(positions.size());
			index.updateUndefined();
			return index;
		}

		BigDecimal getZPart(ConfIndex confIndex, int pos1, int rc1) {

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
					if (pmat.getTuple(new RCTuple(pos1, rc1, pos2, rc2, pos3, rc3).sorted())) {
						return BigDecimal.ZERO;
					}

					BigDecimal triple = blute.get(pos1, rc1, pos2, rc2, pos3, rc3);
					if (triple != null) {
						math.mult(triple);
					}
				}
			}

			return math.get();
		}

		BigDecimal optimizeZ(ConfIndex index, MathTools.Optimizer opt) {

			BigMath energy = bigMath().set(1.0);

			// get the score for each undefined position
			for (int i=0; i<index.numUndefined; i++) {
				int pos1 = index.undefinedPos[i];

				// optimize over possible assignments to pos1
				BigDecimal pos1Energy = opt.initBigDecimal();
				for (int rc1 : rcs.get(pos1)) {

					BigMath rc1Energy = bigMath();
					if (isPruned(index, pos1, rc1)) {

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
								BigDecimal triple = blute.get(pos1, rc1, pos2, rc2, pos3, rc3);
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

										// triple are optional
										BigDecimal triple = blute.get(pos1, rc1, pos2, rc2, pos3, rc3);
										if (triple != null) {
											rc2Energy.mult(triple);
										}
									}

									// triples with undefined positions
									for (int k=0; k<j; k++) {
										int pos3 = index.undefinedPos[k];

										// optimize over rcs
										BigDecimal optrc3Energy = opt.initBigDecimal();
										for (int rc3 : rcs.get(pos3)) {
											BigDecimal triple = blute.get(pos1, rc1, pos2, rc2, pos3, rc3);
											if (triple != null) {
												optrc3Energy = opt.opt(optrc3Energy, triple);
											}
										}

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

		boolean isPruned(ConfIndex index, int pos1, int rc1) {

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
					if (pmat.getTuple(new RCTuple(pos1, rc1, pos2, rc2, pos3, rc3).sorted())) {
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

				// count the RCs by sequence
				Map<String,Integer> counts = new HashMap<>();
				for (int rc : rcs.get(pos)) {
					String resType = positions.get(pos).resConfs.get(rc).template.name;
					counts.compute(resType, (rt, current) -> {
						if (current == null) {
							return 1;
						} else {
							return current + 1;
						}
					});
				}

				int minCount = counts.values().stream().mapToInt(v -> v).min().getAsInt();
				int maxCount = counts.values().stream().mapToInt(v -> v).max().getAsInt();

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

			// pick the next pos to assign
			int pos = index.numDefined;

			for (int rc : rcs.get(pos)) {

				// get the zpart
				BigDecimal zpart;
				if (index.numDefined == 0) {
					zpart = blute.factor;
				} else {
					zpart = getZPart(index, pos, rc);
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
