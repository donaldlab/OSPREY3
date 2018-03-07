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

	// TEMP
	public static void main(String[] args) {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A38", "A39", "A40", "A41")) { // [38,44]
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

			// do DEE
			PruningMatrix pmat = new SimpleDEE.Runner()
				.setSinglesThreshold(100.0)
				.setPairsThreshold(100.0)
				.setShowProgress(true)
				.run(confSpace, emat);

			// get a bunch of conformations and energies
			final int numConfs = 1000;
			final File confDBFile = new File("lute-test.conf.db");
			log("computing energies for %d confs...", numConfs);
			List<ConfSearch.EnergiedConf> econfs = new ConfDB(confSpace, confDBFile).use((confdb) -> {

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
				return confEcalc.calcAllEnergies(
					confs,
					true,
					confdb.new ConfTable("lute-test")
				);
			});

			log("\nLUTE:\n");

			final double maxAllowedResidual = 0.01;

			// compute LUTE matrix for pair tuples
			LUTE lute = new LUTE(confSpace);
			lute.addUnprunedPairTuples(pmat);
			double residual = lute.fit(confEcalc);
			if (residual > maxAllowedResidual) {
				throw new Error(String.format("LUTE residual (%.4f) greater than max allowed (%.4f)", residual, maxAllowedResidual));
			}
			EnergyMatrix lutemat = lute.makeEnergyMatrix();

			// compare conf energies
			for (ConfSearch.EnergiedConf econf : econfs) {
				double luteEnergy = lutemat.confE(econf.getAssignments());
				log("conf %20s   score %9.4f      energy %9.4f   gap %7.4f      LUTE energy %9.4f   diff %7.4f",
					Arrays.toString(econf.getAssignments()),
					econf.getScore(),
					econf.getEnergy(),
					econf.getEnergy() - econf.getScore(),
					luteEnergy,
					luteEnergy - econf.getEnergy()
				);
			}
		}
	}
}
