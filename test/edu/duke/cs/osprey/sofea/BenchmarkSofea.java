package edu.duke.cs.osprey.sofea;


import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.*;

import java.io.File;
import java.math.BigDecimal;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

import static edu.duke.cs.osprey.tools.Log.log;


public class BenchmarkSofea {

	public static void main(String[] args) {

		ForcefieldParams ffparams = new ForcefieldParams();
		File tempDir = new File("/tmp/benchmarkSofea");
		tempDir.mkdirs();

		// use the new templates, cuz why not
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
			.clearTemplateCoords()
			.addTemplateCoords(FileTools.readFile("template coords.v2.txt"))
			.build();

		// define design flexibility [68,73]
		Map<String, List<String>> designFlex = new HashMap<>();
		// unavoidable clash at A68. don't use ARG, or sub something smaller
		//designFlex.put("A68", Arrays.asList(Strand.WildType /* arg=34 */));
		designFlex.put("A69", Arrays.asList(Strand.WildType /* ser=18 *//*, "THR", "LEU", "ILE", "VAL", "ALA", "GLY", "CYS"*/));
		designFlex.put("A70", Arrays.asList(Strand.WildType /* gly=1 *//*, "ALA", "VAL", "LEU", "ILE", "CYS"*/));
		designFlex.put("A71", Arrays.asList(Strand.WildType /* lys=27 */));
		designFlex.put("A72", Arrays.asList(Strand.WildType /* gln=9 */));
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
		Function<List<Strand>, SimpleConfSpace> makeConfSpace = (strands) ->
			new SimpleConfSpace.Builder().addStrands(strands).build();
		MultiStateConfSpace confSpace = new MultiStateConfSpace
			.Builder("complex", makeConfSpace.apply(Arrays.asList(design, target)))
			.build();

		log("seq space: %s", confSpace.seqSpace);

		SofeaLab.ConfEcalcFactory makeConfEcalc = (simpleConfSpace, ecalc, amat) ->
			new ConfEnergyCalculator.Builder(simpleConfSpace, ecalc)
				.setEnergyPartition(EnergyPartition.Traditional)
				.build();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.makeCpu(1)) // single-threaded for benchmarking
			.build()) {

			Sofea sofea = new Sofea.Builder(confSpace)
				.setSweepIncrement(1)
				.setSeqDBFile(new File(tempDir, "sofea.seqdb"))
				.setFringeDBLowerFile(new File(tempDir, "sofea.lower.fringedb"))
				.setFringeDBLowerMiB(32)
				.setFringeDBUpperFile(new File(tempDir, "sofea.upper.fringedb"))
				.setFringeDBUpperMiB(1)
				.configEachState(state -> {

					// always compute emats with all available speed
					EnergyMatrix emat;
					try (EnergyCalculator fastEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
						.setParallelism(Parallelism.makeCpu(4))
						.build()) {

						ConfEnergyCalculator fastConfEcalc = makeConfEcalc.make(state.confSpace, fastEcalc);

						emat = new SimplerEnergyMatrixCalculator.Builder(fastConfEcalc)
							.setCacheFile(new File(tempDir, String.format("sofea.%s.emat", state.name)))
							.setTripleCorrectionThreshold(10.0)
							.build()
							.calcEnergyMatrix();
					}

					ConfEnergyCalculator confEcalc = makeConfEcalc.make(state.confSpace, ecalc);
					return new Sofea.StateConfig(emat, confEcalc, null);
				})
				.build();

			// run pass 1
			sofea.init(true);
			try (FringeDB fringedb = sofea.openFringeDBLower()) {
				try (SeqDB seqdb = sofea.openSeqDB()) {

					Double[] gThresholds = confSpace.states.stream()
						.map(state -> {
							BigExp zSumMax = fringedb.getZSumMax(state);
							if (zSumMax.isNaN()) {
								return null;
							} else {
								return sofea.bcalc.freeEnergyPrecise(zSumMax);
							}
						})
						.toArray(size -> new Double[size]);

					// do an initial step to populate fringedb
					gThresholds[0] += 10.0;
					Stopwatch stopwatch = new Stopwatch().start();
					sofea.pass1(fringedb, seqdb, 0, null, gThresholds, stopwatch, Double.POSITIVE_INFINITY);
					stopwatch.stop();
					log("finished pass 1 step 1 in %s", stopwatch.getTime(2));
					fringedb.finishStep();

					// do another bigger step for benchmarking
					gThresholds[0] += 15.0;
					stopwatch = new Stopwatch().start();
					sofea.pass1(fringedb, seqdb, 0, null, gThresholds, stopwatch, Double.POSITIVE_INFINITY);
					stopwatch.stop();
					log("finished pass 1 step 2 in %s", stopwatch.getTime(2));
				}
			}

			/* TEMP
			MultiStateConfSpace.State state = confSpace.getState("complex");
			Sofea.StateInfo stateInfo = sofea.getStateInfo(state);

			// get the root node bounds
			BigDecimal rootZSumUpper = stateInfo.calcZSumUpper(stateInfo.makeConfIndex(), stateInfo.rcs);
			log("root zSumUpper: %e", rootZSumUpper);

			Stopwatch sw = new Stopwatch().start();
			int numNodes = refineZSumUpper(
				new BigDecimal("2.419422e+0"),
				sofea.zPruneThreshold,
				stateInfo,
				stateInfo.makeConfIndex(),
				rootZSumUpper
			);
			log("finished in %s: num nodes = %d", sw.stop().getTime(2), numNodes);
			*/
		}
	}

	private static int refineZSumUpper(BigDecimal zSumThreshold, BigDecimal zPruneThreshold, Sofea.StateInfo stateInfo, ConfIndex index, BigDecimal zSumUpper) {

		// forget any subtree if it's below the pruning threshold
		if (MathTools.isLessThan(zSumUpper, zPruneThreshold)) {
			return 1;
		}

		// if we're a leaf node, just use the bound
		if (index.isFullyDefined()) {
			return 1;
		}

		// if zSumUpper is too small, add the node to the fringe set
		if (MathTools.isLessThan(zSumUpper, zSumThreshold)) {
			return 1;
		}

		// not a leaf node, recurse
		int numNodes = 1;
		int pos = stateInfo.posPermutation[index.numDefined];
		for (int rc : stateInfo.rcs.get(pos)) {
			index.assignInPlace(pos, rc);
			numNodes += refineZSumUpper(
				zSumThreshold,
				zPruneThreshold,
				stateInfo,
				index,
				stateInfo.calcZSumUpper(index, stateInfo.rcs).toBigDecimal() // TODO: use BigExp further up the stack?
			);
			index.unassignInPlace(pos);
		}

		return numNodes;
	}
}
