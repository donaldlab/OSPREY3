package edu.duke.cs.osprey.sofea;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundFastQueues;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.*;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;

import java.io.File;
import java.math.BigDecimal;
import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.function.Function;

import static edu.duke.cs.osprey.tools.Log.log;


public class SofeaLab {

	public static void main(String[] args) {

		ForcefieldParams ffparams = new ForcefieldParams();
		boolean recalc = false;
		File tempDir = new File("/tmp/sofeaLab");
		tempDir.mkdirs();

		// use the new templates, cuz why not
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
			.clearTemplateCoords()
			.addTemplateCoords(FileTools.readFile("template coords.v2.txt"))
			.build();

		// define design flexibility [68,73]
		Map<String,List<String>> designFlex = new HashMap<>();
		// unavoidable clash at A68. don't use ARG, or sub something smaller
		//designFlex.put("A68", Arrays.asList(Strand.WildType /* arg=34 */));
		//designFlex.put("A69", Arrays.asList(Strand.WildType /* ser=18 *//*, "THR", "LEU", "ILE", "VAL", "ALA", "GLY", "CYS"*/));
		designFlex.put("A70", Arrays.asList(Strand.WildType /* gly=1 *//*, "ALA", "VAL", "LEU", "ILE", "CYS"*/));
		//designFlex.put("A71", Arrays.asList(Strand.WildType /* lys=27 */));
		//designFlex.put("A72", Arrays.asList(Strand.WildType /* gln=9 */));
		designFlex.put("A73", Arrays.asList(Strand.WildType /* leu=5 */));

		// define target flexibility [5,10]
		List<String> targetFlex = Arrays.asList(
			"A5", // lys=27
			"A6", // hie=8
			"A7", // tyr=8
			"A8", // gln=9
			"A9", // phe=4
			"A10" // asn=7
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
		Function<List<Strand>,SimpleConfSpace> makeConfSpace = (strands) ->
			new SimpleConfSpace.Builder().addStrands(strands)
				.setShellDistance(6.0)
				.build();
		MultiStateConfSpace confSpace = new MultiStateConfSpace
			.Builder("complex", makeConfSpace.apply(Arrays.asList(design, target)))
			// TEMP: just complex state for now
			//.addMutableState("design", makeConfSpace.apply(Arrays.asList(design)))
			//.addUnmutableState("target", makeConfSpace.apply(Arrays.asList(target)))
			.build();

		log("seq space: %s", confSpace.seqSpace);

		BiFunction<SimpleConfSpace,EnergyCalculator,ConfEnergyCalculator> makeConfEcalc = (simpleConfSpace, ecalc) ->
			new ConfEnergyCalculator.Builder(simpleConfSpace, ecalc)
				.setEnergyPartition(EnergyPartition.Traditional) // wait for emats is boring...
				//.setEnergyPartition(EnergyPartition.AllOnPairs) // use the tighter lower bounds
				.build();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			//.setParallelism(Parallelism.makeCpu(4))
			.setParallelism(Parallelism.makeCpu(1)) // TEMP: single-threaded for now
			.build()) {

			Sofea sofea = new Sofea.Builder(confSpace)
				.setSweepDivisor(2.0) // TEMP: use a smaller divisor to get more sweep steps
				.configEachState(state -> {

					File ematFile = new File(tempDir, String.format("sofea.%s.emat", state.name));
					File confdbFile = new File(tempDir, String.format("sofea.%s.confdb", state.name));
					if (recalc) {
						ematFile.delete();
						confdbFile.delete();
					}

					ConfEnergyCalculator confEcalc = makeConfEcalc.apply(state.confSpace, ecalc);
					EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
						.setCacheFile(ematFile)
						.build()
						.calcEnergyMatrix();

					return new Sofea.StateConfig(emat, confEcalc, confdbFile);
				})
				.build();

			/* TEMP: brute force the free energies for all sequences
			for (Sequence seq : confSpace.seqSpace.getSequences()) {
				log("%s", seq);
				for (MultiStateConfSpace.State state : confSpace.states) {
					BigDecimal z = sofea.calcZSum(seq, state);
					double g = sofea.bcalc.freeEnergyPrecise(z);
					log("\t%10s   z=%s  g=%.4f", state.name, Log.formatBigLn(z), g);
				}
			}
			*/

			/*
			MultiStateConfSpace.LMFE lmfe = confSpace.lmfe()
				.addPositive("complex")
				.addNegative("design")
				.addNegative("target")
				.build();

			MinLMFE criterion = new MinLMFE(lmfe, 1);

			// do the design!
			sofea.refine(criterion);
			*/


			Consumer<BigDecimalBounds> dumpZ = z -> {
				DoubleBounds g = sofea.bcalc.freeEnergyPrecise(z);
				log("z = %s  d=%9.4f", Log.formatBigLn(z), sofea.bigMath().set(z.upper).sub(z.lower).div(z.upper).get());
				log("g = %s  w=%9.4f", g.toString(4, 9), g.size());
			};


			// compute the free energy of one sequence in the complex state
			Sequence seq = confSpace.seqSpace.makeWildTypeSequence();
			MultiStateConfSpace.State state = confSpace.getState("complex");
			Sofea.StateConfig config = sofea.getConfig(state);

			final double epsilon = 0.14;
			//final double epsilon = 0.0001;

			// using MARK*
			if (false) {
				/* TODO: calc the emat-upper for MARK*
				RCs rcs = seq.makeRCs(state.confSpace);
				MARKStarBoundFastQueues pfunc = new MARKStarBoundFastQueues(
					state.confSpace,
					config.ematUpper,
					config.emat,
					config.confEcalc,
					rcs,
					config.confEcalc.ecalc.parallelism
				);
				pfunc.init(epsilon);
				pfunc.setStabilityThreshold(null);
				pfunc.setReportProgress(true);
				pfunc.setCorrections(new UpdatingEnergyMatrix(state.confSpace, config.emat));
				//pfunc.reduceMinimizations = true or false?
				pfunc.stateName = state.name;

				Stopwatch sw = new Stopwatch().start();
				pfunc.compute();
				log("MARK* pfunc finished in %s", sw.stop().getTime(2));

				dumpZ.accept(new BigDecimalBounds(
					pfunc.getValues().calcLowerBound(),
					pfunc.getValues().calcUpperBound()
				));
				*/
			}

			try (ConfDB confdb = new ConfDB(state.confSpace, new File("pfunc.confdb"))) {
				ConfDB.ConfTable confTable = confdb.new ConfTable("pfunc");

				// using the gradient descent pfunc
				if (false) {
					GradientDescentPfunc pfunc = new GradientDescentPfunc(config.confEcalc);
					RCs rcs = seq.makeRCs(state.confSpace);
					ConfAStarTree astar = new ConfAStarTree.Builder(config.emat, rcs)
						.setTraditional()
						.build();
					pfunc.init(astar, rcs.getNumConformations(), epsilon);
					pfunc.setStabilityThreshold(null);
					pfunc.setReportProgress(true);
					pfunc.setConfTable(confTable);

					Stopwatch sw = new Stopwatch().start();
					pfunc.compute();
					log("GD pfunc finished in %s", sw.stop().getTime(2));
					log("\tnum minimizations: %d", pfunc.getNumConfsEvaluated());

					log("\tZ upper = %s", pfunc.getValues().calcUpperBound());

					dumpZ.accept(new BigDecimalBounds(
						pfunc.getValues().calcLowerBound(),
						pfunc.getValues().calcUpperBound()
					));
				}

				BigDecimal pfuncUpperBound = new BigDecimal("1.489779375480618714133411059887472011543144099753141485887679129E+129");

				// use the confdb to find the "perfect" number of minimizations needed
				if (true) {
					int numConfs = 0;
					BigMath m = sofea.bigMath().set(0);
					for (ConfSearch.EnergiedConf econf : confTable.energiedConfs(ConfDB.SortOrder.Energy)) {

						m.add(sofea.bcalc.calcPrecise(econf.getEnergy()));
						log("\tconf   %4d   [%s]   energy=%9.3f   lower=%9.3f   gap=%9.3f",
							++numConfs,
							Streams.joinToString(state.confSpace.positions, ", ", pos -> {
								int rc = econf.getAssignments()[pos.index];
								return pos.resConfs.get(rc).getRotamerCode();
							}),
							econf.getEnergy(), econf.getScore(), econf.getEnergy() - econf.getScore()
						);

						// are we done yet?
						double delta = sofea.bigMath().set(pfuncUpperBound).sub(m.get()).div(pfuncUpperBound).get().doubleValue();
						if (delta <= epsilon) {
							dumpZ.accept(new BigDecimalBounds(
								m.get(),
								pfuncUpperBound
							));
							log("\tnum minimizations: %d", numConfs);
							break;
						}
					}
				}
			}

			// using SOFEA
			if (true) {
				// clear the confdb
				config.confDBFile.delete();

				MultiStateConfSpace.LMFE lmfe = confSpace.lmfe()
					.addPositive("complex")
					.build();
				double lmfeEnergyWidth = 0.1;
				SequenceLMFE criterion = new SequenceLMFE(seq, lmfe, lmfeEnergyWidth);

				sofea.init(true);
				Stopwatch sw = new Stopwatch().start();
				sofea.refine(criterion);
				log("SOFEA finished in %s", sw.stop().getTime(2));

				try (SeqDB seqdb = sofea.openSeqDB()) {
					dumpZ.accept(seqdb.getSequencedZSumBounds(seq).get(state));
				}
			}
		}
	}
}
