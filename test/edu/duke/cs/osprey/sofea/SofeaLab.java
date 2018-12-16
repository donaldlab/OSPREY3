package edu.duke.cs.osprey.sofea;

import static edu.duke.cs.osprey.tools.Log.logf;
import static edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;

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
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.*;

import java.io.File;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.*;

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
		// unavoidable clash at A68. don't use ARG, or sub something smaller
		//designFlex.put("A68", Arrays.asList(Strand.WildType /* arg */));
		designFlex.put("A69", Arrays.asList(Strand.WildType /* ser */, "THR", "LEU", "ILE", "VAL", "ALA", "GLY", "CYS"));
		designFlex.put("A70", Arrays.asList(Strand.WildType /* gly */, "ALA", "VAL", "LEU", "ILE", "CYS"));
		designFlex.put("A71", Arrays.asList(Strand.WildType /* lys */));
		designFlex.put("A72", Arrays.asList(Strand.WildType /* gln */));
		designFlex.put("A73", Arrays.asList(Strand.WildType /* leu */));

		// define target flexibility [5,10]
		List<String> targetFlex = Arrays.asList(
			"A5", // lys
			"A6", // hie
			"A7", // tyr
			"A8", // gln
			"A9", // phe
			"A10" // asn
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

		// TEMP
		// make a multi-state conf space
		MultiStateConfSpace confSpace = new MultiStateConfSpace
			.Builder("design", new SimpleConfSpace.Builder().addStrands(design).build())
			.addMutableState("complex", new SimpleConfSpace.Builder().addStrands(design, target).build())
			.addUnmutableState("target", new SimpleConfSpace.Builder().addStrands(target).build())
			.build();

		// use the usual affinity optimization objective function
		MinLMFE criterion = new MinLMFE(
			confSpace.lmfe()
				.addPositive("complex")
				.addNegative("design")
				.addNegative("target")
				.build(),
			10,
			new MathContext(16, RoundingMode.HALF_UP)
		);
		//

		/* TEMP: just do the design state
		MultiStateConfSpace confSpace = new MultiStateConfSpace
			.Builder("design", new SimpleConfSpace.Builder().addStrands(design).build())
			.build();

		MinLMFE criterion = new MinLMFE(
			confSpace.lmfe()
				.addPositive("design")
				.build(),
			10,
			new MathContext(16, RoundingMode.HALF_UP)
		);
		*/

		/* TEMP: just do the complex state
		MultiStateConfSpace confSpace = new MultiStateConfSpace
			.Builder("complex", new SimpleConfSpace.Builder().addStrands(design, target).build())
			.build();

		MinLMFE criterion = new MinLMFE(
			confSpace.lmfe()
				.addPositive("complex")
				.build(),
			10,
			new MathContext(16, RoundingMode.HALF_UP)
		);
		*/

		log("seq space: %s", confSpace.seqSpace);

		Sofea sofea;
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.makeCpu(32))
			.build()) {

			// TEMP
			//sofea = new Sofea.Builder(confSpace, null)
			sofea = new Sofea.Builder(confSpace, criterion)
				.setFringeDBMiB(100)
				.configEachState(state -> {

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
						.setThreshold(100.0)
						.setGoldsteinDiffThreshold(50.0)
						.setSinglesPlugThreshold(0.6)
						.setPairsPlugThreshold(0.6)
						//.setTriplesPlugThreshold(0.6)
						.setTransitivePruning(true)
						.setShowProgress(true)
						.run(state.confSpace, emat);

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

					return new Sofea.StateConfig(
						new LUTEConfEnergyCalculator(state.confSpace, luteState),
						pmat
					);
				})
				.make();
		}

		log("\n");

		/* TEMP
		// calc all sequence Z values
		List<Sequence> seqs = new ArrayList<>();
		seqs.add(confSpace.seqSpace.makeWildTypeSequence());
		seqs.addAll(confSpace.seqSpace.getMutants());
		for (Sequence seq : seqs) {
			log("seq %s:", seq);
			for (MultiStateConfSpace.State state : confSpace.states) {
				//Sofea.StateConfig config = sofea.getConfig(state);
				//RCs rcs = new RCs(seq.makeRCs(state.confSpace), config.pmat);
				//BigDecimal z = bruteForcePfuncLuteAStar(config.luteEcalc, rcs);
				BigDecimal z = sofea.calcZ(state, seq);
				log("\t%10s  ln(Z) = %s", state.name, Log.formatBigLn(z));
			}
		}
		*/

		// try SOFEA
		Stopwatch sw = new Stopwatch().start();
		sofea.init(true);
		sofea.refine();
		log("SOFEA:   %9s", sw.stop().getTime(2));

		dump(sofea);
		dumpUnexploredSequences(sofea);
		dumpSequences(sofea);
		try (SeqDB seqdb = sofea.openSeqDB()) {
			criterion.makeResultDoc(seqdb, new File("sofea.md"));
		}
	}

	private static BigDecimal bruteForcePfuncAStar(ConfEnergyCalculator confEcalc, EnergyMatrix emat, RCs rcs) {

		BigMath sum = new BigMath(PartitionFunction.decimalPrecision);
		sum.set(0.0);

		BoltzmannCalculator bcalc = new BoltzmannCalculator(sum.context);

		ConfAStarTree astar = new ConfAStarTree.Builder(emat, rcs)
			.setTraditional()
			.build();
		while (true) {

			ConfSearch.ScoredConf conf = astar.nextConf();
			if (conf == null) {
				break;
			}

			confEcalc.calcEnergyAsync(conf, (epmol) -> {
				sum.add(bcalc.calcPrecise(epmol.getEnergy()));
			});
		}

		return sum.get();
	}

	private static BigDecimal bruteForcePfuncLuteAStar(LUTEConfEnergyCalculator luteEcalc, RCs rcs) {

		BigMath sum = new BigMath(PartitionFunction.decimalPrecision);
		sum.set(0.0);

		BoltzmannCalculator bcalc = new BoltzmannCalculator(sum.context);

		ConfAStarTree astar = new ConfAStarTree.Builder(null, rcs)
			.setLUTE(luteEcalc)
			.build();
		for (ConfSearch.ScoredConf conf : astar.nextConfs(Double.POSITIVE_INFINITY)) {
			sum.add(bcalc.calcPrecise(conf.getScore()));
		}

		return sum.get();
	}

	private static void dump(Sofea sofea) {

		try (SeqDB seqdb = sofea.openSeqDB()) {

			log("Seq DB:");
			log("\tunsequenced state bounds:");
			for (MultiStateConfSpace.State state : seqdb.confSpace.unsequencedStates) {
				BigDecimalBounds z = seqdb.getUnsequencedBound(state);
				Log.log("\t\t%10s=%s", state.name, Log.formatBigLn(z));
			}
			log("\tsequenced state sums:");
			for (Map.Entry<Sequence,SeqDB.SeqInfo> entry : seqdb.getSequencedSums()) {
				Sequence seq = entry.getKey();
				SeqDB.SeqInfo seqInfo = entry.getValue();

				logf("\t\t");
				for (MultiStateConfSpace.State state : seqdb.confSpace.sequencedStates) {
					BigDecimalBounds bounds = seqInfo.z[state.sequencedIndex];
					logf("%10s=%s   ", state.name, Log.formatBigLn(bounds));
				}
				log("[%s]", seq);
			}
		}
	}

	private static void dumpUnexploredSequences(Sofea sofea) {

		try (SeqDB seqdb = sofea.openSeqDB()) {

			log("Seq DB: unexplored sequences:");

			for (Map.Entry<Sequence,SeqDB.SeqInfo> entry : seqdb.getSequencedBounds()) {

				Sequence seq = entry.getKey();
				SeqDB.SeqInfo seqInfo = entry.getValue();

				if (seq.isFullyAssigned()) {
					continue;
				}

				logf("\t");
				for (MultiStateConfSpace.State state : seqdb.confSpace.sequencedStates) {
					logf("%10s=%s   ", state.name, Log.formatBigLn(seqInfo.z[state.sequencedIndex]));

				}
				log("[%s]", seq);
			}
		}
	}

	private static void dumpSequences(Sofea sofea) {

		try (SeqDB seqdb = sofea.openSeqDB()) {

			log("Seq DB: sequences:");
			for (Map.Entry<Sequence,SeqDB.SeqInfo> entry : seqdb.getSequencedBounds()) {

				Sequence seq = entry.getKey();
				SeqDB.SeqInfo seqInfo = entry.getValue();

				if (!seq.isFullyAssigned()) {
					continue;
				}

				logf("\t");
				for (MultiStateConfSpace.State state : seqdb.confSpace.sequencedStates) {
					BigDecimalBounds bounds = seqInfo.z[state.sequencedIndex];
					logf("%10s=%s d=%.4f   ", state.name, Log.formatBigLn(bounds), bounds.delta(sofea.mathContext));
				}
				log("[%s]", seq);
			}
		}
	}
}
