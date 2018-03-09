package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static edu.duke.cs.osprey.tools.Log.log;

public class LUTEPlayground {

	public static void main(String[] args) {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		//List<String> resNums = Arrays.asList("A38", "A39", "A40"); // this works great on pairs
		List<String> resNums = Arrays.asList("A38", "A39", "A40", "A41"); // NOTE: need triples to get a good LUTE fit here
		for (String resNum : resNums) { // [38,44]
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType, "VAL", "LEU")
				.addWildTypeRotamers()
				.setContinuous();
		}

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.make(5, 0, 0))
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

			// compute conventional energy matrix
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(new File("lute-test.emat.dat"))
				.build()
				.calcEnergyMatrix();

			// prune clashes
			PruningMatrix pmat = new SimpleDEE.Runner()
				.setSinglesThreshold(100.0)
				.setPairsThreshold(100.0)
				.setShowProgress(true)
				.run(confSpace, emat);

			// get a bunch of conformations and energies
			final int numConfs = 10;
			final File confDBFile = new File("lute-test.conf.db");
			log("computing energies for %d confs...", numConfs);
			new ConfDB(confSpace, confDBFile).use((confdb) -> {
				ConfDB.ConfTable confTable = confdb.new ConfTable("lute-test");

				ConfSearch astar = new ConfAStarTree.Builder(emat, pmat)
					.setTraditional()
					.build();

				// gather the confs
				List<ConfSearch.ScoredConf> confs = new ArrayList<>();
				for (int i=0; i<numConfs; i++) {
					ConfSearch.ScoredConf conf = astar.nextConf();
					if (conf == null) {
						break;
					}
					confs.add(conf);
				}

				// compute the energies
				List<ConfSearch.EnergiedConf> econfs = confEcalc.calcAllEnergies(confs, true, confTable);

				log("\nLUTE:\n");

				final int minSamplesPerTuple = 10;

				// compute LUTE matrix for pair tuples
				LUTE lute = new LUTE(confSpace);
				lute.addUnprunedPairTuples(pmat);
				LUTE.Errors errors = lute.fit(confEcalc, confTable, minSamplesPerTuple);

				// compare conf energies
				for (ConfSearch.EnergiedConf econf : econfs) {
					double luteEnergy = lute.emat.confE(econf.getAssignments());
					log("conf %20s   score %9.4f      energy %9.4f   gap %7.4f      LUTE energy %9.4f   diff %7.4f",
						Arrays.toString(econf.getAssignments()),
						econf.getScore(),
						econf.getEnergy(),
						econf.getEnergy() - econf.getScore(),
						luteEnergy,
						luteEnergy - econf.getEnergy()
					);
				}
			});
		}
	}
}
