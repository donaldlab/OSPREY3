package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Log;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Supplier;

import static edu.duke.cs.osprey.tools.Log.log;

public class LUTEPlayground {

	public static void main(String[] args) {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		List<String> mutableResNums = Arrays.asList("A38", "A39");
		//List<String> mutableResNums = Arrays.asList("A38", "A39", "A40"); // this works great on pairs
		//List<String> mutableResNums = Arrays.asList("A38", "A39", "A40", "A41"); // NOTE: need triples to get a good LUTE fit here
		for (String resNum : mutableResNums) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType, "VAL", "LEU")
				.addWildTypeRotamers()
				.setContinuous();
		}

		// add some flexibility
		List<String> flexibleResNums = Arrays.asList("A40", "A41", "A42", "A43", "A44");
		for (String resNum : flexibleResNums) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.make(8, 0, 0))
			//.setParallelism(Parallelism.make(8, 3, 16))
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

			// compute conventional energy matrix
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(new File("lute-test.emat.dat"))
				.build()
				.calcEnergyMatrix();

			// run DEE (super important for good LUTE fits!!)
			PruningMatrix pmat = new SimpleDEE.Runner()
				.setSinglesThreshold(100.0)
				.setPairsThreshold(100.0)
				.setGoldsteinDiffThreshold(20.0)
				.setShowProgress(true)
				.run(confSpace, emat);
			//PruningMatrix pmat = new PruningMatrix(confSpace);

			Supplier<ConfSearch> astarFactory = () -> new ConfAStarTree.Builder(emat, pmat)
				.setTraditional()
				.build();

			LUTE lute;

			final File confDBFile = new File("lute-test.conf.db");
			try (ConfDB confdb = new ConfDB(confSpace, confDBFile)) {
				ConfDB.ConfTable confTable = confdb.new ConfTable("lute-test");

				// compute energies for some low-energy confs
				final int numConfs = 10;
				log("computing energies for %d confs...", numConfs);
				List<ConfSearch.ScoredConf> confs = new ArrayList<>();
				{
					ConfSearch astar = astarFactory.get();
					for (int i=0; i<numConfs; i++) {
						ConfSearch.ScoredConf conf = astar.nextConf();
						if (conf == null) {
							break;
						}
						confs.add(conf);
					}
				}
				List<ConfSearch.EnergiedConf> econfs = confEcalc.calcAllEnergies(confs, true, confTable);

				log("\nLUTE:\n");

				final int randomSeed = 12345;
				final int minSamplesPerTuple = 10;

				// compute LUTE matrix for pair tuples
				lute = new LUTE(confSpace);
				lute.addUnprunedPairTuples(pmat);
				//lute.addUnprunedTripleTuples(pmat);
				LUTE.Errors errors = lute.fit(confEcalc, confTable, minSamplesPerTuple, randomSeed);

				// attempt to count the conf space by enumeration
				int confSpaceSize = 0;
				boolean isExhausted = false;
				{
					ConfSearch astar = astarFactory.get();
					for (confSpaceSize=0; confSpaceSize<errors.getNumSamples()*2; confSpaceSize++) {
						if (astar.nextConf() == null) {
							isExhausted = true;
							break;
						}
					}
				}
				if (isExhausted) {
					log("conf space (after all pruning) has EXACTLY %d confs", confSpaceSize);
					if (confSpaceSize == errors.getNumSamples()) {
						log("LUTE sampling has entirely exhausted the conf space, can't sample any more confs for training or testing");
					}
				} else {
					log("conf space (after all pruning) has at least %d confs", confSpaceSize);
				}

				// compare conf energies
				for (ConfSearch.EnergiedConf econf : econfs) {
					double luteEnergy = lute.emat.confE(econf.getAssignments());
					log("conf %30s   score %9.4f      energy %9.4f   gap %7.4f      LUTE energy %9.4f   diff %7.4f",
						Arrays.toString(econf.getAssignments()),
						econf.getScore(),
						econf.getEnergy(),
						econf.getEnergy() - econf.getScore(),
						luteEnergy,
						luteEnergy - econf.getEnergy()
					);
				}
			}

			final double pfuncEpsilon = 0.1;

			/* NOPE

			// compute a partition function the old-fashioned way
			log("\nplain old pfunc:\n");
			Stopwatch oldSw = new Stopwatch().start();
			PartitionFunction oldPfunc = new SimplePartitionFunction(astarFactory.get(), confEcalc);
			oldPfunc.init(pfuncEpsilon);
			//oldPfunc.setReportProgress(true);
			oldPfunc.compute();
			log("done in %s   %s", oldSw.stop().getTime(2), oldPfunc.makeResult());

			// use LUTE to do the same thing
			log("\nLUTE-powered pfunc:\n");
			Stopwatch luteSw = new Stopwatch().start();
			LUTEConfEnergyCalculator luteEcalc = new LUTEConfEnergyCalculator(confSpace, ecalc, lute.emat);
			PartitionFunction lutePfunc = new SimplePartitionFunction(astarFactory.get(), luteEcalc);
			lutePfunc.init(pfuncEpsilon);
			//oldPfunc.setReportProgress(true);
			lutePfunc.compute();
			log("done in %s   %s", luteSw.stop().getTime(2), lutePfunc.makeResult());
			*/
		}
	}
}
