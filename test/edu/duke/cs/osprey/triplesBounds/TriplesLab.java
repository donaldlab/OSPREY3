package edu.duke.cs.osprey.triplesBounds;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.order.SequentialAStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;
import java.util.Arrays;
import java.util.function.Function;

import static edu.duke.cs.osprey.tools.Log.log;

public class TriplesLab {

	private static interface AStarFactory {
		ConfAStarTree make();
	}

	public static void main(String[] args) {

		// get an arbitrary conf space
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9")) {
			strand.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		}
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		RCs rcs = new RCs(confSpace);

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			 .setParallelism(Parallelism.makeCpu(8))
			 .build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				//.setEnergyPartition(EnergyPartition.Traditional)
				.setEnergyPartition(EnergyPartition.AllOnPairs)
				.build();

			AStarFactory astarFactory;
			if (false) {

				// compute an energy matrix
				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
					.setCacheFile(new File("temp/emat.dat"))
					.build()
					.calcEnergyMatrix();

				astarFactory = () -> new ConfAStarTree.Builder(emat, rcs)
					.setTraditional()
					.build();

			} else {

				// compute triples energy
				TriplesEnergy triplesEnergy = new TriplesEnergy(confSpace);
				triplesEnergy.calculateOrCache(ecalc, new File("temp/triples.dat"));
				assert (triplesEnergy.isFullyDefined());

				astarFactory = () -> new ConfAStarTree.Builder(null, rcs)
					.setCustom(
						//new DynamicHMeanAStarOrder(),
						//new StaticScoreHMeanAStarOrder(),
						new SequentialAStarOrder(),
						triplesEnergy.new GScorer(),
						triplesEnergy.new HScorer()
					)
					.build();
			}

			// analyze bounds looseness
			analyzeConfs(astarFactory, confEcalc, (econf) -> econf.getScore());

			// compuate a (single-threaded, confdb-backed) partition function
			try (ConfDB confDB = new ConfDB(confSpace, new File("temp/pfunc.confdb"))) {
				ConfDB.ConfTable confTable = confDB.new ConfTable("foo");

				try (EnergyCalculator ecalcSingle = new EnergyCalculator.Builder(confSpace, new ForcefieldParams()).build()) {
					ConfEnergyCalculator confEcalcSingle = new ConfEnergyCalculator.Builder(confSpace, ecalcSingle).build();

					final double epsilon = 0.5;
					GradientDescentPfunc pfunc = new GradientDescentPfunc(confEcalcSingle);
					pfunc.init(astarFactory.make(), rcs.getNumConformations(), epsilon);
					pfunc.setReportProgress(true);
					pfunc.setConfTable(confTable);
					pfunc.compute();

					log("%d conformations evaluated", pfunc.getNumConfsEvaluated());
				}
			}
		}
	}

	private static void analyzeConfs(AStarFactory astarFactory, ConfEnergyCalculator confEcalc, Function<ConfSearch.EnergiedConf,Double> scorer) {

		ConfAStarTree astar = astarFactory.make();

		final int numConfs = 10;
		double rmse = 0.0;

		for (int i=0; i<numConfs; i++) {

			ConfSearch.ScoredConf conf = astar.nextConf();
			if (conf == null) {
				break;
			}

			ConfSearch.EnergiedConf econf = confEcalc.calcEnergy(conf);

			double score = scorer.apply(econf);
			double error = econf.getEnergy() - score;
			log("conf %d %s    score %.6f    energy %.6f    error %.6f", i, Conf.toString(econf.getAssignments()), score, econf.getEnergy(), error);
			error = Math.abs(error);
			rmse += error*error;
		}

		rmse = Math.sqrt(rmse/numConfs);
		log("top %d confs RMSE: %.6f", numConfs, rmse);
	}
}
