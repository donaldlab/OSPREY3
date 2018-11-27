package edu.duke.cs.osprey.sofea;

import static edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import static edu.duke.cs.osprey.tools.MathTools.BigIntegerBounds;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.NewalgLab;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.lute.*;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.*;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;


public class SofeaLab {

	public static void main(String[] args) {

		ForcefieldParams ffparams = new ForcefieldParams();

		// use the new templates, cuz why not
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
			.clearTemplateCoords()
			.addTemplateCoords(FileTools.readFile("template coords.v2.txt"))
			.build();

		boolean recalc = false;

		// define design flexibility [68,73]
		Map<String,List<String>> designFlex = new HashMap<>();
		//designFlex.put("A68", Arrays.asList(Strand.WildType /* arg */));
		designFlex.put("A69", Arrays.asList(Strand.WildType /* ser */, "THR"));
		designFlex.put("A70", Arrays.asList(Strand.WildType /* gly */, "ALA"));
		//designFlex.put("A71", Arrays.asList(Strand.WildType /* lys */));
		//designFlex.put("A72", Arrays.asList(Strand.WildType /* gln */));
		designFlex.put("A73", Arrays.asList(Strand.WildType /* leu */));

		// define target flexibility [5,10]
		List<String> targetFlex = Arrays.asList(
			"A5", // lys
			"A6" // hie
			// "A7" // tyr
			// "A8" // gln
			// "A9" // phe
			// "A10" // asn
		);

		// build strands
		Molecule pdb = PDBIO.readResource("/1CC8.ss.pdb");
		Strand design = new Strand.Builder(pdb)
			.setTemplateLibrary(templateLib)
			.setResidues("A68", "A73")
			.build();
		for (Map.Entry<String,List<String>> entry : designFlex.entrySet()) {
			design.flexibility.get(entry.getKey())
				.setLibraryRotamers(entry.getValue())
				.addWildTypeRotamers()
				.setContinuous();
		}
		Strand target = new Strand.Builder(pdb)
			.setTemplateLibrary(templateLib)
			.setResidues("A2", "A67")
			.build();
		for (String resNum : targetFlex) {
			target.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}

		// make a multi-state conf space
		MultiStateConfSpace confSpace = new MultiStateConfSpace
			.Builder("design", new SimpleConfSpace.Builder().addStrands(design).build())
			.addSequencedState("complex", new SimpleConfSpace.Builder().addStrands(design, target).build())
			.addState("target", new SimpleConfSpace.Builder().addStrands(target).build())
			.build();

		log("seq space: %s", confSpace.seqSpace);

		List<Sofea.StateConfig> stateConfigs = new ArrayList<>();
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			for (MultiStateConfSpace.State state : confSpace.states) {

				File ematFile = new File(String.format("sofea.%s.emat", state.name));
				File pmatFile = new File(String.format("sofea.%s.pmat", state.name));
				File luteFile = new File(String.format("sofea.%s.lute", state.name));
				if (recalc) {
					ematFile.delete();
					pmatFile.delete();
					luteFile.delete();
				}

				ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(state.confSpace, ecalc).build();
				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
					.setCacheFile(ematFile)
					.build()
					.calcEnergyMatrix();
				PruningMatrix pmat = new SimpleDEE.Runner()
					.setCacheFile(pmatFile)
					.setParallelism(ecalc.parallelism)
					.setThreshold(null)
					.setSinglesPlugThreshold(0.6)
					.setPairsPlugThreshold(0.6)
					.setTriplesPlugThreshold(0.6)
					.setShowProgress(true)
					.run(state.confSpace, null);

				// do LUTE stuff
				LUTEState luteState;
				if (luteFile.exists()) {
					luteState = LUTEIO.read(luteFile);
					log("read LUTE state from file: %s", luteFile.getAbsolutePath());
				} else {
					try (ConfDB confdb = new ConfDB(state.confSpace)) {
						ConfDB.ConfTable confTable = confdb.new ConfTable("lute");

						final int randomSeed = 12345;
						final LUTE.Fitter fitter = LUTE.Fitter.OLSCG;
						final double maxOverfittingScore = 1.5;
						final double maxRMSE = 0.1;

						// compute LUTE fit
						LUTE lute = new LUTE(state.confSpace);
						ConfSampler sampler = new RandomizedDFSConfSampler(state.confSpace, pmat, randomSeed);
						lute.sampleTuplesAndFit(confEcalc, emat, pmat, confTable, sampler, fitter, maxOverfittingScore, maxRMSE);
						lute.reportConfSpaceSize(pmat);

						luteState = new LUTEState(lute.getTrainingSystem());
						LUTEIO.write(luteState, luteFile);
						log("wrote LUTE state to file: %s", luteFile.getAbsolutePath());
					}
				}

				stateConfigs.add(new Sofea.StateConfig(
					new LUTEConfEnergyCalculator(state.confSpace, luteState),
					pmat
				));
			}

			log("\n");

			// calc all sequence Z values
			{
				List<Sequence> seqs = new ArrayList<>();
				seqs.add(confSpace.seqSpace.makeWildTypeSequence());
				seqs.addAll(confSpace.seqSpace.getMutants());
				for (Sequence seq : seqs) {
					log("seq %s:", seq);
					for (MultiStateConfSpace.State state : confSpace.states) {
						Sofea.StateConfig config = stateConfigs.get(state.index);
						NewalgLab.PfuncCalc pcalc = new NewalgLab.PfuncCalc(config.luteEcalc, config.pmat, seq.makeRCs(state.confSpace));
						log("\t%10s  ln(Z) = %s", state.name, dump(pcalc.calc()));
					}
				}
			}

			// try SOFEA
			{
				Sofea sofea = new Sofea(confSpace, stateConfigs);
				Stopwatch sw = new Stopwatch().start();
				sofea.design();
				log("SOFEA:   %9s", sw.stop().getTime(2));
			}
		}
	}

	public static String dump(ConfIndex index) {
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

	public static String dump(BigDecimal val) {
		if (val == null) {
			return "null";
		}
		return String.format("%9.4f", new BoltzmannCalculator(PartitionFunction.decimalPrecision).ln(val));
	}

	public static String dump(MathTools.BigDecimalBounds values) {
		return String.format("[%-9s,%9s]",
			dump(values.lower),
			dump(values.upper)
		);
	}

	public static double getBoundsDelta(BigDecimalBounds bounds) {
		return new BigMath(PartitionFunction.decimalPrecision)
			.set(bounds.upper)
			.sub(bounds.lower)
			.div(bounds.upper)
			.get()
			.doubleValue();
	}

	public static class Sofea {

		public static class StateConfig {

			public final LUTEConfEnergyCalculator luteEcalc;
			public final PruningMatrix pmat;

			public StateConfig(LUTEConfEnergyCalculator luteEcalc, PruningMatrix pmat) {
				this.luteEcalc = luteEcalc;
				this.pmat = pmat;
			}
		}

		public final MultiStateConfSpace confSpace;
		public final List<StateConfig> stateConfigs;

		// NOTE: don't need much precision for most math, but need lots of precision for seqdb math
		// TODO: externalize these as config options
		private final MathContext mathContext = new MathContext(32, RoundingMode.HALF_UP);
		private final MathContext seqdbMathContext = new MathContext(256, RoundingMode.HALF_UP);

		private final List<StateInfo> stateInfos;

		public Sofea(MultiStateConfSpace confSpace, List<StateConfig> stateConfigs) {

			this.confSpace = confSpace;
			this.stateConfigs = stateConfigs;

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

		private BigMath bigMath() {
			return new BigMath(mathContext);
		}

		public void design() {

			try (SeqDB seqdb = new SeqDB(confSpace, seqdbMathContext)) {

				FringeDB fringedb = new FringeDB(confSpace);

				// process the root node for each state
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
					{
						ConfIndex index = stateInfo.makeConfIndex();
						fringedb.add(state, index, rootBound, stateInfo.blute.factor);
						seqdb.addSeqZ(state, index, rootBound);
					}
				}

				// TEMP
				log("seed");
				log("\tfringe size: %d", fringedb.writeSize());

				final double zmaxFactor = Math.pow(Math.E, 4.0);

				// TEMP: try a few operations
				for (int i=0; i<10; i++) {

					// reduce zmax each iteration
					List<BigDecimal> zmax = confSpace.states.stream()
						.map(state -> bigMath()
							.set(fringedb.getZMax(state))
							.div(zmaxFactor)
							.get()
						)
						.collect(Collectors.toList());

					// TEMP
					log("op %3d", i);
					for (MultiStateConfSpace.State state : confSpace.states) {
						 log("\tzmax=%s  /f=%s", dump(fringedb.getZMax(state)), dump(zmax.get(state.index)));
					}

					long[] numProcessed = { 0 };
					long oldFringeSize = fringedb.writeSize();

					fringedb.sweep((state, index, bounds, zpath) -> {

						seqdb.subSeqZ(state, index, bounds);

						boolean wasProcessed = design(state, index, bounds, zmax.get(state.index), zpath, fringedb, seqdb);

						if (wasProcessed) {
							numProcessed[0]++;
						}
					});

					// TEMP
					log("\tprocessed %d/%d, fringe size: %d", numProcessed[0], oldFringeSize, fringedb.writeSize());

					// stop if we ran out of fringe
					if (fringedb.writeSize() <= 0) {
						break;
					}
				}

				// TEMP
				seqdb.dump();
				seqdb.dumpPartialSequences();
				seqdb.dumpSequences();
			}
		}

		private boolean design(MultiStateConfSpace.State state, ConfIndex index, BigDecimalBounds bounds, BigDecimal zmax, BigDecimal zpath, FringeDB fringedb, SeqDB seqdb) {

			// skip this tree if it's too small
			if (MathTools.isLessThan(bounds.upper, zmax)) {
				fringedb.add(state, index, bounds, zpath);
				seqdb.addSeqZ(state, index, bounds);
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
					// TODO: can optimize seq creation
					index.assignInPlace(pos, rc);
					seqdb.addSeqZ(state, index, zpathrc);
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
					design(state, index, boundsrc, zmax, zpathrc, fringedb, seqdb);
					index.unassignInPlace(pos);
				}
			}

			return true;
		}

		private class StateInfo {

			final MultiStateConfSpace.State state;
			final NewalgLab.BoltzmannLute blute;
			final PruningMatrix pmat;
			final RCs rcs;
			final List<SimpleConfSpace.Position> positions;

			StateInfo(MultiStateConfSpace.State state) {
				this.state = state;
				StateConfig config = stateConfigs.get(state.index);
				this.blute = new NewalgLab.BoltzmannLute(config.luteEcalc, mathContext);
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

			BigIntegerBounds count(ConfIndex index) {

				BigIntegerBounds count = new BigIntegerBounds(BigInteger.ONE, BigInteger.ONE);

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
		}
	}
}
