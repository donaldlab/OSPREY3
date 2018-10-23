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
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static edu.duke.cs.osprey.tools.Log.log;


public class SofeaLab {

	public static void main(String[] args) {

		ForcefieldParams ffparams = new ForcefieldParams();

		// use the new templates, cuz why not
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
			.clearTemplateCoords()
			.addTemplateCoords(FileTools.readFile("template coords.v2.txt"))
			.build();

		// define flexibility
		Map<String,List<String>> flex = new HashMap<>();
		flex.put("A2", Arrays.asList(Strand.WildType /* ala */, "GLY"));
		flex.put("A3", Arrays.asList(Strand.WildType /* glu */, "ASP"));
		flex.put("A4", Arrays.asList(Strand.WildType /* ile */, "LEU"));
		flex.put("A5", Arrays.asList(Strand.WildType /* lys */, "ARG"));
		flex.put("A6", Arrays.asList(Strand.WildType /* hie */, "PHE"));
		flex.put("A7", Arrays.asList(Strand.WildType /* tyr */, "PHE"));
		boolean recalc = false;

		// TODO: try to analyze the partial order of seq leaves and trees?
		// TODO: how good of a partial order is good enough to stop refining?

		// make a conf space
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb"))
			.setTemplateLibrary(templateLib)
			.build();
		for (Map.Entry<String,List<String>> entry : flex.entrySet()) {
			strand.flexibility.get(entry.getKey())
				.setLibraryRotamers(entry.getValue())
				.addWildTypeRotamers()
				.setContinuous();
		}
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		log("seq space: %s", confSpace.seqSpace);

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
				.setSinglesPlugThreshold(0.6)
				.setPairsPlugThreshold(0.6)
				.setTriplesPlugThreshold(0.6)
				.setShowProgress(true)
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

			log("\n");

			/* TEMP
			// calc all sequence Z values
			{
				List<Sequence> seqs = new ArrayList<>();
				seqs.add(confSpace.makeWildTypeSequence());
				seqs.addAll(confSpace.seqSpace.getMutants());
				for (Sequence seq : seqs) {
					PfuncCalc pcalc = new PfuncCalc(luteEcalc, pmat, seq.makeRCs(confSpace));
					log("seq %s:  ln(Z) = %s", seq, dump(pcalc.calc()));
				}
			}
			*/

			// calc multi-sequence Z tight bounds using the new alg
			{
				NewalgLab.PfuncCalc pcalc = new NewalgLab.PfuncCalc(luteEcalc, pmat);
				Stopwatch sw = new Stopwatch().start();
				MathTools.BigDecimalBounds Zbounds = pcalc.calcMultiSequence();
				log("NewAlg:      %9s   ln(Z) in %s", sw.stop().getTime(2), dump(Zbounds));
			}

			/* TEMP
			// estimate multi-sequence Z bounds using the new alg
			//for (double pruneFactor : Arrays.asList(1.0, 0.9, 0.5, 0.1, 0.01, 0.0)) {
			{ double pruneFactor = 0.5;
				int depth = confSpace.positions.size();
				PfuncCalc pcalc = new PfuncCalc(luteEcalc, pmat);
				Stopwatch sw = new Stopwatch().start();
				MathTools.BigDecimalBounds Zbounds = pcalc.estimateMultiSequence(depth, pruneFactor);
				log("NewAlg %4.2f:   %9s   ln(Z) in %s", pruneFactor, sw.stop().getTime(2), dump(Zbounds));
			}
			*/

			// try SOFEA
			{
				Sofea sofea = new Sofea(luteEcalc, pmat);
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

		public final LUTEConfEnergyCalculator luteEcalc;
		public final PruningMatrix pmat;

		private final MathContext mathContext = PartitionFunction.decimalPrecision;

		private final NewalgLab.BoltzmannLute blute;
		private final int numPos;
		private final RCs rcs;

		public Sofea(LUTEConfEnergyCalculator luteEcalc, PruningMatrix pmat) {
			this(luteEcalc, pmat, new RCs(luteEcalc.confSpace));
		}

		public Sofea(LUTEConfEnergyCalculator luteEcalc, PruningMatrix pmat, RCs rcs) {

			this.luteEcalc = luteEcalc;
			this.pmat = pmat;
			this.rcs = rcs;

			this.blute = new NewalgLab.BoltzmannLute(luteEcalc, mathContext);
			this.numPos = luteEcalc.confSpace.positions.size();
		}

		private BigMath bigMath() {
			return new BigMath(mathContext);
		}

		public void design() {

			if (luteEcalc.confSpace.positions.size() == 1) {

				// TODO
				throw new UnsupportedOperationException("TODO");

			} else {

				// get a multi-sequence Z bound on the root node
				MathTools.BigDecimalBounds rootBound;
				{
					ConfIndex index = new ConfIndex(luteEcalc.confSpace.positions.size());
					index.updateUndefined();

					BigDecimal minZ = optimizeZ(index, MathTools.Optimizer.Minimize);
					BigDecimal maxZ = optimizeZ(index, MathTools.Optimizer.Maximize);
					MathTools.BigIntegerBounds size = count(index);

					rootBound = new MathTools.BigDecimalBounds(
						bigMath()
							.set(minZ)
							.mult(size.lower)
							.mult(blute.factor)
							.get(),
						bigMath()
							.set(maxZ)
							.mult(size.upper)
							.mult(blute.factor)
							.get()
					);
				}

				try (SeqDB seqdb = new SeqDB(luteEcalc.confSpace, mathContext)) {

					// init the fringe with the root node
					FringeDB fringedb = new FringeDB(luteEcalc.confSpace);
					{
						ConfIndex index = new ConfIndex(luteEcalc.confSpace.positions.size());
						index.updateUndefined();
						fringedb.add(index, rootBound, blute.factor);
						seqdb.addSeqZ(index, rootBound);
					}

					// TEMP: try a few operations
					for (int i=0; i<20; i++) {

						// reduce zmax each iteration
						BigDecimal zmax = bigMath()
							.set(fringedb.getZMax())
							.div(Math.E)
							.get();
						log("op %3d  zmax=%s  /2=%s root=%s", i, dump(fringedb.getZMax()), dump(zmax), dump(rootBound));

						final BigDecimal fzmax = zmax; // silly Java closures...
						fringedb.sweep((index, bounds, zpath) -> {

							seqdb.subSeqZ(index, bounds);

							design(index, bounds, fzmax, zpath, fringedb, seqdb);
						});

						// TEMP
						log("fringe size: %d", fringedb.writeSize());
					}

					// TEMP
					seqdb.dumpPartialSequences();
					seqdb.dumpSequences();
				}
			}
		}

		private void design(ConfIndex index, BigDecimalBounds bounds, BigDecimal zmax, BigDecimal zpath, FringeDB fringedb, SeqDB seqdb) {

			// skip this tree if it's too small
			if (MathTools.isLessThan(bounds.upper, zmax)) {
				fringedb.add(index, bounds, zpath);
				seqdb.addSeqZ(index, bounds);
				return;
			}

			// otherwise, recurse

			boolean isRoot = index.numDefined == 0;
			boolean isRCLeaf = index.numDefined + 1 == index.numPos;

			int pos = index.numDefined;
			for (int rc : rcs.get(pos)) {

				// update the zpath with this RC
				BigDecimal zpathrc = zpath;
				if (!isRoot) {

					BigDecimal zrc = getZPart(index, pos, rc);

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
					seqdb.addSeqZ(index, zpathrc);
					index.unassignInPlace(pos);

				} else {

					// get the subtree bounds
					index.assignInPlace(pos, rc);
					BigDecimal minZ = optimizeZ(index, MathTools.Optimizer.Minimize);
					BigDecimal maxZ = optimizeZ(index, MathTools.Optimizer.Maximize);
					MathTools.BigIntegerBounds count = count(index);
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
					design(index, boundsrc, zmax, zpathrc, fringedb, seqdb);
					index.unassignInPlace(pos);
				}
			}
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

		private BigDecimal optimizeZ(ConfIndex index, MathTools.Optimizer opt) {

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

		private BigIntegerBounds count(ConfIndex index) {

			BigIntegerBounds count = new BigIntegerBounds(BigInteger.ONE, BigInteger.ONE);

			for (int i=0; i<index.numUndefined; i++) {
				int pos = index.undefinedPos[i];

				// count the RCs by sequence
				Map<String,Integer> counts = new HashMap<>();
				for (int rc : rcs.get(pos)) {
					String resType = luteEcalc.confSpace.positions.get(pos).resConfs.get(rc).template.name;
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
