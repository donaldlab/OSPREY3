package edu.duke.cs.osprey.kstar;


import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.lute.*;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.*;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;


public class NewalgLab {

	public static void main(String[] args) {

		ForcefieldParams ffparams = new ForcefieldParams();

		// use the new templates, cuz why not
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
			.clearTemplateCoords()
			.addTemplateCoords(FileTools.readFile("template coords.v2.txt"))
			.build();

		// make a conf space
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb"))
			.setTemplateLibrary(templateLib)
			.build();
		List<String> resNums = Arrays.asList(
			"A2",
			"A3",
			"A4",
			"A5",
			"A6",
			"A7"
		);
		for (String resNum : resNums) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		boolean recalc = false;
		File ematFile = new File("newalg.emat");
		File pmatFile = new File("newalg.pmat");
		File luteFile = new File("newalg.lute");
		if (recalc) {
			ematFile.delete();
			pmatFile.delete();
			luteFile.delete();
		}

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(ematFile)
				.build()
				.calcEnergyMatrix();
			PruningMatrix pmat = new SimpleDEE.Runner()
				.setCacheFile(pmatFile)
				.setParallelism(ecalc.parallelism)
				.setThreshold(null)
				//.setSinglesGoldsteinDiffThreshold(20.0)
				//.setPairsGoldsteinDiffThreshold(20.0)
				//.setTriplesGoldsteinDiffThreshold(20.0)
				.setSinglesPlugThreshold(0.6)
				.setPairsPlugThreshold(0.6)
				.setTriplesPlugThreshold(0.6)
				.setShowProgress(true)
				//.run(confSpace, emat);
				.run(confSpace, null);

			// do LUTE stuff
			LUTEState luteState;
			if (luteFile.exists()) {
				luteState = LUTEIO.read(luteFile);
				log("read LUTE state from file: %s", luteFile.getAbsolutePath());
			} else {
				try (ConfDB confdb = new ConfDB(confSpace)) {
					ConfDB.ConfTable confTable = confdb.new ConfTable("lute");

					final int randomSeed = 12345;
					final LUTE.Fitter fitter = LUTE.Fitter.OLSCG;
					final double maxOverfittingScore = 1.5;
					final double maxRMSE = 0.1;

					// compute LUTE fit
					LUTE lute = new LUTE(confSpace);
					ConfSampler sampler = new RandomizedDFSConfSampler(confSpace, pmat, randomSeed);
					lute.sampleTuplesAndFit(confEcalc, emat, pmat, confTable, sampler, fitter, maxOverfittingScore, maxRMSE);
					lute.reportConfSpaceSize(pmat);

					luteState = new LUTEState(lute.getTrainingSystem());
					LUTEIO.write(luteState, luteFile);
					log("wrote LUTE state to file: %s", luteFile.getAbsolutePath());
				}
			}
			LUTEConfEnergyCalculator luteEcalc = new LUTEConfEnergyCalculator(confSpace, luteState);

			// TODO: prune all the LUTE zero-energy tuples?

			log("\n");

			// calc the LUTE approximation to Z exactly (ie, as well as LUTE could ever do)
			{
				Stopwatch sw = new Stopwatch().start();
				BigDecimal Z = bruteForcePfuncLuteAStar(luteEcalc, pmat);
				log("LUTE A*:    %9s   ln(Z) = %s", sw.stop().getTime(2), dump(Z));
			}

			// estimate Z using the LUTE pfunc caluclator
			{
				Stopwatch sw = new Stopwatch().start();
				PartitionFunction.Values Z = lutePfunc(luteEcalc, pmat, 0.9);
				log("LUTE Pcalc: %9s   ln(Z) in %s", sw.stop().getTime(2), dump(Z));
			}

			// calc Z using the new alg
			{
				PfuncCalc pcalc = new PfuncCalc(luteEcalc, pmat);
				Stopwatch sw = new Stopwatch().start();
				BigDecimal Z = pcalc.calc();
				log("NewAlg:     %9s   ln(Z) = %s", sw.stop().getTime(2), dump(Z));
			}

			// estimate Z using the new alg
			double pruneFactor = 0.0;
			for (int depth=0; depth<=confSpace.positions.size(); depth++) {
			//for (int depth=0; depth<=3; depth++) {
			//{ int depth = 4;
				PfuncCalc pcalc = new PfuncCalc(luteEcalc, pmat);
				Stopwatch sw = new Stopwatch().start();
				PartitionFunction.Values Z = pcalc.estimate(depth, pruneFactor);
				log("NewAlg %d:   %9s   log(Z) in %s", depth, sw.stop().getTime(2), dump(Z));
			}

			// TEMP
			//new PfuncCalc(luteEcalc, pmat).test();
		}
	}

	private static String dump(ConfIndex index) {
		StringBuilder buf = new StringBuilder();
		buf.append('[');
		for (int i=0; i<index.numDefined; i++) {
			if (i > 0) {
				buf.append(", ");
			}
			buf.append(index.definedPos[i]);
			buf.append('=');
			buf.append(index.definedRCs[i]);
		}
		buf.append(']');
		return buf.toString();
	}

	private static String dump(BigDecimal val) {
		if (val == null) {
			return "null";
		}
		// TEMP
		val = new BigMath(PartitionFunction.decimalPrecision)
			.set(val)
			//.add(BigDecimal.ONE)
			.get();
		return String.format("%7.2f", new BoltzmannCalculator(PartitionFunction.decimalPrecision).ln(val));
		//return String.format("%e", val);
	}

	private static String dump(PartitionFunction.Values values) {
		// TEMP
		if (false) {
			return dumpFreeEnergy(values);
		}
		return String.format("[%-9s,%9s] d=%.6f",
			dump(values.calcLowerBound()),
			dump(values.calcUpperBound()),
			values.getEffectiveEpsilon()
		);
	}

	private static String dumpFreeEnergy(PartitionFunction.Values values) {
		BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
		double lo = bcalc.freeEnergyPrecise(values.calcUpperBound());
		double hi = bcalc.freeEnergyPrecise(values.calcLowerBound());
		double delta = (lo - hi)/hi;
		return String.format("[%-7.2f,%7.2f] d=%.6f", lo, hi, delta);
	}

	private static BigDecimal bruteForcePfuncAStar(ConfEnergyCalculator confEcalc, EnergyMatrix emat) {

		BigDecimal sum = BigDecimal.ZERO;
		BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

		ConfAStarTree astar = new ConfAStarTree.Builder(emat, confEcalc.confSpace)
			.setTraditional()
			.build();
		List<ConfSearch.ScoredConf> confs = astar.nextConfs(Double.POSITIVE_INFINITY);
		List<ConfSearch.EnergiedConf> econfs = confEcalc.calcAllEnergies(confs);
		for (ConfSearch.EnergiedConf conf : econfs) {
			sum = sum.add(bcalc.calc(conf.getEnergy()));
		}

		return sum;
	}

	private static BigDecimal bruteForcePfuncLuteAStar(LUTEConfEnergyCalculator luteEcalc, PruningMatrix pmat) {

		BigDecimal sum = BigDecimal.ZERO;
		BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

		ConfAStarTree astar = new ConfAStarTree.Builder(null, pmat)
			.setLUTE(luteEcalc)
			.build();
		for (ConfSearch.ScoredConf conf : astar.nextConfs(Double.POSITIVE_INFINITY)) {
			sum = sum.add(bcalc.calc(conf.getScore()));
		}

		return sum;
	}

	private static PartitionFunction.Values lutePfunc(LUTEConfEnergyCalculator luteEcalc, PruningMatrix pmat, double epsilon) {

		ConfAStarTree astar = new ConfAStarTree.Builder(null, pmat)
			.setLUTE(luteEcalc)
			.build();

		LUTEPfunc pfunc = new LUTEPfunc(luteEcalc);
		pfunc.init(astar, new RCs(luteEcalc.confSpace).getNumConformations(), epsilon);
		pfunc.setStabilityThreshold(null);
		//pfunc.setReportProgress(true);
		pfunc.compute();

		return pfunc.getValues();
	}

	public static class BoltzmannLute {

		public final LUTEConfEnergyCalculator luteEcalc;

		private final BigDecimal factor;
		private final BigDecimal[] values;

		public BoltzmannLute(LUTEConfEnergyCalculator luteEcalc, MathContext mathContext) {

			this.luteEcalc = luteEcalc;
			this.values = new BigDecimal[luteEcalc.tuples.size()];

			// pre-calculate all the boltzmann-weighted values
			BoltzmannCalculator bcalc = new BoltzmannCalculator(mathContext);
			this.factor = bcalc.calcPrecise(luteEcalc.state.tupleEnergyOffset);
			for (int t=0; t<luteEcalc.tuples.size(); t++) {
				this.values[t] = bcalc.calcPrecise(luteEcalc.state.tupleEnergies[t]);
				/* TEMP
				log("\tLUTE[%4d] = e=%16.12f %s z=%e  tuple=%s",
					t,
					luteEcalc.state.tupleEnergies[t],
					luteEcalc.state.tupleEnergies[t] == 0.0 ? "Z" : " ",
					values[t],
					luteEcalc.tuples.get(t)
				);
				*/
			}

			// TEMP
			//log("offset: %.4f   %s", luteEcalc.state.tupleEnergyOffset, dump(factor));
		}

		public boolean has(int pos, int rc) {
			return luteEcalc.hasTuple(pos, rc);
		}

		public boolean has(int pos1, int rc1, int pos2, int rc2) {
			return luteEcalc.hasTuple(pos1, rc1, pos2, rc2);
		}

		public BigDecimal get(int pos, int rc) {
			return get(luteEcalc.tuples.getIndex(pos, rc));
		}

		public BigDecimal get(int pos1, int rc1, int pos2, int rc2) {
			return get(luteEcalc.tuples.getIndex(pos1, rc1, pos2, rc2));
		}

		public BigDecimal get(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {
			return get(luteEcalc.tuples.getIndex(pos1, rc1, pos2, rc2, pos3, rc3));
		}

		private BigDecimal get(Integer t) {
			if (t == null) {
				return null;
			} else {
				return values[t];
			}
		}
	}

	/* TODO: LUTE sort of normalizes the tuple energies,
		can we do the log sum of exponentials with doubles?
		see: https://jblevins.org/log/log-sum-exp
	 */
	public static class PfuncCalc {

		public final LUTEConfEnergyCalculator luteEcalc;
		public final PruningMatrix pmat;

		private final MathContext mathContext = PartitionFunction.decimalPrecision;

		private final BoltzmannLute blute;
		private final int numPos;
		private final RCs rcs;

		public PfuncCalc(LUTEConfEnergyCalculator luteEcalc, PruningMatrix pmat) {

			this.luteEcalc = luteEcalc;
			this.pmat = pmat;

			this.blute = new BoltzmannLute(luteEcalc, mathContext);
			this.numPos = luteEcalc.confSpace.positions.size();
			this.rcs = new RCs(luteEcalc.confSpace);
		}

		private BigMath bigMath() {
			return new BigMath(mathContext);
		}

		public BigDecimal calc() {

			if (luteEcalc.confSpace.positions.size() == 1) {

				// for tiny conf spaces, assume only all singles present
				BigMath math = bigMath()
					.set(0.0);
				for (int rc : rcs.get(0)) {
					BigDecimal val = blute.get(0, rc);
					if (val != null) {
						math.add(val);
					}
				}
				return math
					.mult(blute.factor)
					.get();

			} else {

				// otherwise, assume no singles, all pairs, and maybe some triples present
				ConfIndex index = new ConfIndex(luteEcalc.confSpace.positions.size());
				index.updateUndefined();
				return calc(index).multiply(blute.factor);
			}
		}

		public PartitionFunction.Values estimate(int maxDepth, double pruneFactor) {

			PartitionFunction.Values values = PartitionFunction.Values.makeFullRange();

			if (luteEcalc.confSpace.positions.size() == 1) {

				// for tiny confspaces, just calc the whole pfunc
				values.qstar = calc();
				values.qprime = BigDecimal.ZERO;

			} else {

				ConfIndex index = new ConfIndex(luteEcalc.confSpace.positions.size());
				index.updateUndefined();

				if (maxDepth == 0) {

					BigDecimal leaf = findLeaf(index);
					BigDecimal maxZ = optimizeZ(index, MathTools.Optimizer.Maximize);
					BigInteger size = count(index);

					values.qstar = bigMath()
						.set(leaf)
						.mult(blute.factor)
						.get();
					values.qprime = bigMath()
						.set(maxZ)
						.mult(size)
						.mult(blute.factor)
						.sub(values.qstar)
						.get();

				} else {

					// get the bound width at the root node
					BigDecimal rootBoundWidth = bigMath()
						.set(optimizeZ(index, MathTools.Optimizer.Maximize))
						.mult(count(index))
						.sub(findLeaf(index))
						.get();

					// scale the bound width by the pruning factor
					// but do the scaling in log space, ie exp(ln(w + 1)*f) - 1
					BoltzmannCalculator bcalc = new BoltzmannCalculator(mathContext);
					BigDecimal maxBoundWidth = bigMath()
						.set(
							bcalc.exp(
								bcalc.ln(
									bigMath()
										.set(rootBoundWidth)
										.add(1.0)
										.get()
								)*pruneFactor
							)
						)
						.sub(1.0)
						.get();

					maxDepth = Math.min(maxDepth, index.numPos);
					MathTools.BigDecimalBounds Z = estimate(index, maxDepth, maxBoundWidth);

					values.qstar = bigMath()
						.set(Z.lower)
						.mult(blute.factor)
						.get();
					values.qprime = bigMath()
						.set(Z.upper)
						.mult(blute.factor)
						.sub(values.qstar)
						.get();
				}
			}

			return values;
		}

		private MathTools.BigDecimalBounds estimate(ConfIndex index, int maxDepth, BigDecimal maxBoundWidth) {

			BigMath lomath = bigMath().set(0.0);
			BigMath himath = bigMath().set(0.0);

			boolean isRoot = index.numDefined == 0;
			boolean isRCLeaf = index.numDefined + 1 == index.numPos;

			int pos = index.numDefined;
			for (int rc : rcs.get(pos)) {

				// get the Z part added by this rc, if possible
				BigDecimal zpart;
				if (isRoot) {

					// I AM ROOT
					zpart = BigDecimal.ONE;

				} else {

					zpart = getZPart(index, pos, rc);

					// this subtree contributes nothing to Z
					if (MathTools.isZero(zpart)) {
						continue;
					}
				}

				if (isRCLeaf) {

					assert (!isRoot);

					// hit bottom,  just add up the Z part
					lomath.add(zpart);
					himath.add(zpart);

				} else {

					// get the subtree sum bounds
					index.assignInPlace(pos, rc);
					BigDecimal leaf = findLeaf(index);
					BigDecimal maxZ = optimizeZ(index, MathTools.Optimizer.Maximize);
					BigInteger count = count(index);
					index.unassignInPlace(pos);

					// how wide is the bound?
					BigDecimal boundWidth = bigMath()
						.set(maxZ)
						.mult(count)
						.sub(leaf)
						.get();
					boolean isBoundTooWide = MathTools.isGreaterThanOrEqual(boundWidth, maxBoundWidth);

					// TODO: check for short circuits here?

					// should we recurse here?
					// TODO: check bound width
					boolean canGoDeeper = index.numDefined < maxDepth - 1;
					if (canGoDeeper && isBoundTooWide) {

						// yup, recurse
						index.assignInPlace(pos, rc);
						MathTools.BigDecimalBounds subbounds = estimate(index, maxDepth, maxBoundWidth);
						index.unassignInPlace(pos);

						lomath.add(bigMath()
							.set(zpart)
							.mult(subbounds.lower)
							.get()
						);
						himath.add(bigMath()
							.set(zpart)
							.mult(subbounds.upper)
							.get()
						);

					} else {

						// nope, just use the bound
						lomath.add(bigMath()
							.set(zpart)
							.mult(leaf)
							.get()
						);
						himath.add(bigMath()
							.set(zpart)
							.mult(maxZ)
							.mult(count)
							.get()
						);
					}
				}
			}

			return new MathTools.BigDecimalBounds(lomath.get(), himath.get());
		}

		private BigInteger count(ConfIndex index) {
			BigInteger count = BigInteger.ONE;
			for (int i=0; i<index.numUndefined; i++) {
				int pos = index.undefinedPos[i];
				count = count.multiply(BigInteger.valueOf(rcs.get(pos).length));
			}
			return count;
		}

		private BigDecimal optimizeZ(ConfIndex index, MathTools.Optimizer opt) {

			// TEMP
			//log("%s Z", opt);

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

								// TEMP
								//log("\t\tundef rc2 %d,%d = %s", pos2, rc2, dump(rc2Energy.get()));

								optrc2Energy = opt.opt(optrc2Energy, rc2Energy.get());
							}

							// TEMP
							//log("\t\tpos2 %d = %s", pos2, dump(optrc2Energy));

							rc1Energy.mult(optrc2Energy);
						}
					}

					// TEMP
					//log("\tundef rc1 %d,%d = %s", pos1, rc1, dump(rc1Energy.get()));

					pos1Energy = opt.opt(pos1Energy, rc1Energy.get());
				}

				// TEMP
				//log("\tpos1 %d = %s", pos1, dump(pos1Energy));

				assert (MathTools.isFinite(pos1Energy));

				energy.mult(pos1Energy);
			}

			return energy.get();
		}

		private boolean isPruned(ConfIndex index, int pos1, int rc1) {

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

		private static class RCInfo implements Comparable<RCInfo> {

			public final int rc;
			public final BigDecimal zpart;

			public RCInfo(int rc, BigDecimal zpart) {
				this.rc = rc;
				this.zpart = zpart;
			}

			@Override
			public int compareTo(RCInfo other) {
				// negate so RCs are sorted in weakly descending order
				return -this.zpart.compareTo(other.zpart);
			}
		}

		private BigDecimal findLeaf(ConfIndex index) {

			boolean isRoot = index.numDefined == 0;
			boolean isLeaf = index.numDefined == index.numPos;
			boolean recurse = index.numDefined + 1 < index.numPos;

			// doesn't work on leaf nodes
			assert (!isLeaf);

			// pick the next pos to assign
			int pos = index.numDefined;

			List<RCInfo> rcInfos;
			if (isRoot) {
				// no Z parts for RC at root node, so don't bother sorting
				rcInfos = Arrays.stream(this.rcs.get(pos))
					.mapToObj(rc -> new RCInfo(rc, BigDecimal.ONE))
					.collect(Collectors.toList());
			} else {
				// sort RCs by tuple Z part, to try to bias search towards biggest Z
				rcInfos = Arrays.stream(this.rcs.get(pos))
					.mapToObj(rc -> new RCInfo(rc, getZPart(index, pos, rc)))
					.filter(info -> MathTools.isPositive(info.zpart))
					.sorted()
					.collect(Collectors.toList());
			}

			for (RCInfo rcInfo : rcInfos) {
				int rc = rcInfo.rc;

				BigMath rcmath = bigMath().set(rcInfo.zpart);

				if (recurse) {

					// have positions left to assign, so recurse
					index.assignInPlace(pos, rc);
					BigDecimal subval = findLeaf(index);
					index.unassignInPlace(pos);

					// short circuit for efficiency
					if (MathTools.isZero(subval)) {
						continue;
					}

					rcmath.mult(subval);
				}

				return rcmath.get();
			}

			return BigDecimal.ZERO;
		}

		private BigDecimal calc(ConfIndex index) {

			// cannot work on leaf nodes
			assert (index.numDefined < index.numPos);

			BigMath math = bigMath().set(0.0);

			// pick the next pos to assign
			int pos = index.numDefined;

			for (int rc : rcs.get(pos)) {

				if (index.numDefined == 0) {

					// first pos: no singles, so just recurse
					index.assignInPlace(pos, rc);
					math.add(calc(index));
					index.unassignInPlace(pos);

				} else {

					// not first pos: have pairs, so get additional boltzmann-weighed energy
					BigDecimal e = getZPart(index, pos, rc);

					// short circuit for efficiency
					if (MathTools.isZero(e)) {
						continue;
					}

					if (index.numDefined < numPos - 1) {

						// have positions left to assign, so recurse
						index.assignInPlace(pos, rc);
						BigDecimal subval = calc(index);
						index.unassignInPlace(pos);

						math.add(bigMath()
							.set(subval)
							.mult(e)
							.get()
						);

					} else {

						// all positions assigned, just add the additional sum
						math.add(e);
					}
				}
			}

			assert (math.get() != null);
			return math.get();
		}

		private BigDecimal getZPart(ConfIndex confIndex, int pos1, int rc1) {

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
	}
}
