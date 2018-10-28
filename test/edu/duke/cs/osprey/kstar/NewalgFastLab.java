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
import java.util.*;
import java.util.function.Function;

import static edu.duke.cs.osprey.tools.Log.log;


public class NewalgFastLab {

	public static void main(String[] args) {

		ForcefieldParams ffparams = new ForcefieldParams();

		// use the new templates, cuz why not
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
			.clearTemplateCoords()
			.addTemplateCoords(FileTools.readFile("template coords.v2.txt"))
			.build();

		// define flexibility
		List<String> flex = Arrays.asList(
			"A2",
			"A3",
			"A4",
			"A5",
			"A6",
			"A7",
			"A8",
			"A9"
		);
		boolean recalc = false;

		// make a conf space
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb"))
			.setTemplateLibrary(templateLib)
			.build();
		for (String resNum : flex) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

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
				.setGoldsteinDiffThreshold(20.0)
				//.setPlugThreshold(0.6)
				.setShowProgress(true)
				.run(confSpace, emat);

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
				NewalgLab.PfuncCalc pcalc = new NewalgLab.PfuncCalc(luteEcalc, pmat);
				Stopwatch sw = new Stopwatch().start();
				BigDecimal Z = pcalc.calc();
				log("NewAlg:     %9s   ln(Z) = %s", sw.stop().getTime(2), dump(Z));
			}

			// calc Z using the new fast alg
			{
				PfuncCalcFast pcalc = new PfuncCalcFast(luteEcalc, pmat);
				Stopwatch sw = new Stopwatch().start();
				double Z = pcalc.calc();
				log("NewFast:    %9s   ln(Z) = %14.8f", sw.stop().getTime(2), Math.log(Z));
			}
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

	public static String dump(BigDecimal val) {
		if (val == null) {
			return "null";
		}
		return String.format("%14.8f", new BoltzmannCalculator(PartitionFunction.decimalPrecision).ln(val));
	}

	private static String dumpScaled(MathTools.BigDecimalBounds values, BigDecimal scale) {
		return String.format("[%-9s,%9s]",
			dump(new BigMath(PartitionFunction.decimalPrecision)
				.set(values.lower)
				.mult(scale)
				.get()
			),
			dump(new BigMath(PartitionFunction.decimalPrecision)
				.set(values.upper)
				.mult(scale)
				.get()
			)
		);
	}

	private static String dump(PartitionFunction.Values values) {
		return String.format("[%-9s,%9s] d=%.6f",
			dump(values.calcLowerBound()),
			dump(values.calcUpperBound()),
			values.getEffectiveEpsilon()
		);
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

	public static class BoltzmannLuteFast {

		public final LUTEConfEnergyCalculator luteEcalc;

		public final double factor;
		private final double[] values;

		public BoltzmannLuteFast(LUTEConfEnergyCalculator luteEcalc) {

			this.luteEcalc = luteEcalc;
			this.values = new double[luteEcalc.tuples.size()];

			// E_lute = factor*(sum E_tuple)
			// c + ln(exp(e1 - c) + exp(e2 - c))

			Function<Double,Double> exp = (e) ->
				Math.exp(-e/BoltzmannCalculator.constRT);

			// pre-calculate all the boltzmann-weighted values
			this.factor = exp.apply(luteEcalc.state.tupleEnergyOffset);
			for (int t=0; t<luteEcalc.tuples.size(); t++) {
				this.values[t] = exp.apply(luteEcalc.state.tupleEnergies[t]);
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
			//log("offset: %.4f   %s", luteEcalc.state.tupleEnergyOffset, factor);
		}

		public boolean has(int pos, int rc) {
			return luteEcalc.hasTuple(pos, rc);
		}

		public boolean has(int pos1, int rc1, int pos2, int rc2) {
			return luteEcalc.hasTuple(pos1, rc1, pos2, rc2);
		}

		public Double get(int pos, int rc) {
			return get(luteEcalc.tuples.getIndex(pos, rc));
		}

		public Double get(int pos1, int rc1, int pos2, int rc2) {
			return get(luteEcalc.tuples.getIndex(pos1, rc1, pos2, rc2));
		}

		public Double get(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {
			return get(luteEcalc.tuples.getIndex(pos1, rc1, pos2, rc2, pos3, rc3));
		}

		private Double get(Integer t) {
			if (t == null) {
				return null;
			} else {
				return values[t];
			}
		}
	}

	public static class PfuncCalcFast {

		public final LUTEConfEnergyCalculator luteEcalc;
		public final PruningMatrix pmat;

		private final BoltzmannLuteFast blute;
		private final int numPos;
		private final RCs rcs;

		public PfuncCalcFast(LUTEConfEnergyCalculator luteEcalc, PruningMatrix pmat) {

			this.luteEcalc = luteEcalc;
			this.pmat = pmat;
			this.rcs = new RCs(luteEcalc.confSpace);

			this.blute = new BoltzmannLuteFast(luteEcalc);
			this.numPos = luteEcalc.confSpace.positions.size();
		}

		public double calc() {

			ConfIndex index = new ConfIndex(luteEcalc.confSpace.positions.size());
			index.updateUndefined();
			return calc(index) * blute.factor;
		}

		private double calc(ConfIndex index) {

			// cannot work on leaf nodes
			assert (index.numDefined < index.numPos);

			double out = 0.0;

			// pick the next pos to assign
			int pos = index.numDefined;
			for (int rc : rcs.get(pos)) {

				if (index.numDefined == 0) {

					// first pos: no singles, so just recurse
					index.assignInPlace(pos, rc);
					out += calc(index);
					index.unassignInPlace(pos);

				} else {

					// not first pos: have pairs, so get additional boltzmann-weighed energy
					double zpart = getZPart(index, pos, rc);

					// short circuit for efficiency
					if (zpart == 0.0) {
						continue;
					}

					if (index.numDefined < numPos - 1) {

						// have positions left to assign, so recurse
						index.assignInPlace(pos, rc);
						out += calc(index)*zpart;
						index.unassignInPlace(pos);

					} else {

						// all positions assigned, just add the additional sum
						out += zpart;
					}
				}
			}

			return out;
		}

		private double getZPart(ConfIndex confIndex, int pos1, int rc1) {

			assert (confIndex.numDefined > 0);

			if (pmat.getOneBody(pos1, rc1)) {
				return 0.0;
			}

			double out = 1.0;
			for (int i=0; i<confIndex.numDefined; i++) {
				int pos2 = confIndex.definedPos[i];
				int rc2 = confIndex.definedRCs[i];

				if (pmat.getPairwise(pos1, rc1, pos2, rc2)) {
					return 0.0;
				}

				out *= blute.get(pos1, rc1, pos2, rc2);

				for (int j=0; j<i; j++) {
					int pos3 = confIndex.definedPos[j];
					int rc3 = confIndex.definedRCs[j];

					if (pmat.getPairwise(pos1, rc1, pos3, rc3)) {
						return 0.0;
					}
					if (pmat.getPairwise(pos2, rc2, pos3, rc3)) {
						return 0.0;
					}
					if (pmat.getTuple(new RCTuple(pos1, rc1, pos2, rc2, pos3, rc3).sorted())) {
						return 0.0;
					}

					Double triple = blute.get(pos1, rc1, pos2, rc2, pos3, rc3);
					if (triple != null) {
						out *= triple;
					}
				}
			}

			return out;
		}
	}
}
