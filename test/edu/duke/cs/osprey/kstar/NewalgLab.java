package edu.duke.cs.osprey.kstar;


import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
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
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.List;

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

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.build()
				.calcEnergyMatrix();
			PruningMatrix pmat = new SimpleDEE.Runner()
				.setParallelism(ecalc.parallelism)
				.setThreshold(null)
				//.setSinglesGoldsteinDiffThreshold(20.0)
				//.setPairsGoldsteinDiffThreshold(20.0)
				//.setTriplesGoldsteinDiffThreshold(20.0)
				.setSinglesPlugThreshold(0.4)
				.setPairsPlugThreshold(0.4)
				.setTriplesPlugThreshold(0.4)
				.setShowProgress(true)
				//.run(confSpace, emat);
				.run(confSpace, null);

			// fit a LUTE model
			LUTEState luteState;
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
			}
			LUTEConfEnergyCalculator luteEcalc = new LUTEConfEnergyCalculator(confSpace, luteState);

			log("\n");

			// calc the LUTE approximation to Z exactly (ie, as well as LUTE could ever do)
			{
				Stopwatch sw = new Stopwatch().start();
				BigDecimal Z = bruteForcePfuncLuteAStar(luteEcalc, pmat);
				log("LUTE A*:    %9s   log(Z) = %s", sw.stop().getTime(2), KStarScore.scoreToLog10String(Z));
			}

			// estimate Z using the LUTE pfunc caluclator
			{
				Stopwatch sw = new Stopwatch().start();
				PartitionFunction.Values Z = lutePfunc(luteEcalc, pmat, 0.1);
				log("LUTE Pcalc: %9s   log(Z) in %s", sw.stop().getTime(2), dump(Z));
			}

			// estimate Z using the new alg
			{
				PfuncCalc pcalc = new PfuncCalc(luteEcalc, pmat);
				Stopwatch sw = new Stopwatch().start();
				PartitionFunction.Values Z = pcalc.calc();
				log("NewAlg:     %9s   log(Z) in %s", sw.stop().getTime(2), dump(Z));
			}
		}
	}

	private static String dump(PartitionFunction.Values values) {
		return String.format("[%-9s,%9s] d=%.3f",
			KStarScore.scoreToLog10String(values.calcLowerBound()),
			KStarScore.scoreToLog10String(values.calcUpperBound()),
			values.getEffectiveEpsilon()
		);
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
			this.factor = bcalc.calc(luteEcalc.state.tupleEnergyOffset);
			for (int t=0; t<luteEcalc.tuples.size(); t++) {
				this.values[t] = bcalc.calcPrecise(luteEcalc.state.tupleEnergies[t]);
			}
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

			//hscorer = new LUTEHScorer(luteEcalc);
		}

		private BigMath bigMath() {
			return new BigMath(mathContext);
		}

		public PartitionFunction.Values calc() {

			PartitionFunction.Values values = PartitionFunction.Values.makeFullRange();

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
				values.qstar = math
					.mult(blute.factor)
					.get();

				// estimates for tiny conf spaces are perfect
				values.qprime = BigDecimal.ZERO;

			} else {

				// otherwise, assume no singles, all pairs, and maybe some triples present

				// TEMP
				values.qstar = helper(new RCTuple()).multiply(blute.factor);
				values.qprime = BigDecimal.ZERO;
			}

			return values;
		}

		private BigDecimal helper(RCTuple tuple) {

			BigMath math = bigMath().set(0.0);

			// pick the next pos to assign
			int pos = tuple.size();

			for (int rc : rcs.get(pos)) {

				if (tuple.size() == 0) {

					// TEMP
					//log("pos=%d, rc=%d", pos, rc);

					// first pos, no singles, so just recurse
					tuple.pos.add(pos);
					tuple.RCs.add(rc);
					math.add(helper(tuple));
					tuple.pos.remove(tuple.pos.size() - 1);
					tuple.RCs.remove(tuple.RCs.size() - 1);

				} else {

					// not first pos, have pairs, so get additional bit
					BigDecimal additionalSum = getAdditionalSum(tuple, pos, rc);
					// TEMP
					//log("tuple=%s   pos=%d, rc=%d   sum=%e", tuple, pos, rc, additionalSum);
					if (additionalSum == null) {
						// rc not compatible with this tuple, skip it
						continue;
					}

					if (tuple.size() < numPos - 1) {

						// have positions left to assign, so recurse
						tuple.pos.add(pos);
						tuple.RCs.add(rc);
						math.add(math.group()
							.set(helper(tuple))
							.mult(additionalSum)
							.get()
						);
						tuple.pos.remove(tuple.pos.size() - 1);
						tuple.RCs.remove(tuple.RCs.size() - 1);

					} else {

						// all positions assigned, just add the additional sum
						math.add(additionalSum);
					}
				}
			}

			return math.get();
		}

		private BigDecimal getAdditionalSum(RCTuple tuple, int pos1, int rc1) {

			assert (tuple.size() > 0);

			if (pmat.getOneBody(pos1, rc1)) {
				return null;
			}

			BigMath math = bigMath().set(1.0);
			for (int i=0; i<tuple.size(); i++) {
				int pos2 = tuple.pos.get(i);
				int rc2 = tuple.RCs.get(i);

				if (pmat.getPairwise(pos1, rc1, pos2, rc2)) {
					return null;
				}

				math.mult(blute.get(pos1, rc1, pos2, rc2));

				for (int j=0; j<i; j++) {
					int pos3 = tuple.pos.get(j);
					int rc3 = tuple.RCs.get(j);

					if (pmat.getPairwise(pos1, rc1, pos3, rc3)) {
						return null;
					}
					if (pmat.getPairwise(pos2, rc2, pos3, rc3)) {
						return null;
					}
					if (pmat.getTuple(new RCTuple(pos1, rc1, pos2, rc2, pos3, rc3).sorted())) {
						return null;
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
