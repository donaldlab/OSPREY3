package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Supplier;

import static edu.duke.cs.osprey.tools.Log.log;


public class LUTELab {

	public static void main(String[] args) {
		gmecMain();
		//kstarMain();
	}

	private static SimpleConfSpace make1CC8ConfSpace() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		List<String> mutableResNums = Arrays.asList("A38", "A39");
		for (String resNum : mutableResNums) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType, "VAL", "LEU")
				.addWildTypeRotamers()
				.setContinuous();
		}

		// add some flexibility
		List<String> flexibleResNums = Arrays.asList("A40", "A41", "A42", "A43", "A44", "A45");
		for (String resNum : flexibleResNums) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}

		return new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();
	}

	public static void kstarMain() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make1GUA11();

		train("protein", confSpaces.protein);
		train("ligand", confSpaces.ligand);
		train("complex", confSpaces.complex);
	}

	private static void train(String name, SimpleConfSpace confSpace) {

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(8))
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

			// compute energy matrix
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(new File(String.format("LUTE.%s.emat.dat", name)))
				.build()
				.calcEnergyMatrix();

			// run DEE (super important for good LUTE fits!!)
			PruningMatrix pmat = new SimpleDEE.Runner()
				.setSinglesThreshold(100.0)
				.setPairsThreshold(100.0)
				.setGoldsteinDiffThreshold(5.0)
				.setShowProgress(true)
				.run(confSpace, emat);

			final File confDBFile = new File(String.format("LUTE.%s.conf.db", name));
			try (ConfDB confdb = new ConfDB(confSpace, confDBFile)) {
				ConfDB.ConfTable confTable = confdb.new ConfTable("lute");

				log("\nLUTE:\n");

				final int randomSeed = 12345;
				//final LUTE.Fitter fitter = LUTE.Fitter.LASSO;
				final LUTE.Fitter fitter = LUTE.Fitter.OLSCG;
				final double maxOverfittingScore = 1.5;
				final double maxRMSE = 0.1;

				confEcalc.resetCounters();

				// compute LUTE fit
				LUTE lute = new LUTE(confSpace);
				ConfSampler sampler = new UniformConfSampler(confSpace, pmat, randomSeed);
				lute.sampleTuplesAndFit(confEcalc, pmat, confTable, sampler, fitter, maxOverfittingScore, maxRMSE);
				lute.reportConfSpaceSize(pmat);
				lute.save(new File(String.format("LUTE.%s.dat", name)));
			}
		}
	}

	public static void gmecMain() {

		SimpleConfSpace confSpace = make1CC8ConfSpace();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.make(8, 0, 0))
			//.setParallelism(Parallelism.make(8, 3, 16))
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

			// compute conventional energy matrix
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(new File("LUTE.emat.dat"))
				.build()
				.calcEnergyMatrix();

			// run DEE (super important for good LUTE fits!!)
			PruningMatrix pmat = new SimpleDEE.Runner()
				.setSinglesThreshold(100.0)
				.setPairsThreshold(100.0)
				.setGoldsteinDiffThreshold(10.0)
				.setShowProgress(true)
				.run(confSpace, emat);
				//PruningMatrix pmat = new PruningMatrix(confSpace);

			File file = new File("LUTE.dat");

			train(confSpace, ecalc, confEcalc, emat, pmat, file);
			//astar(confSpace, ecalc, confEcalc, emat, pmat, file);
			//gmec(confSpace, ecalc, confEcalc, emat, pmat, file);
			//pfunc(confSpace, ecalc, confEcalc, emat, pmat, file);
		}
	}

	private static void train(SimpleConfSpace confSpace, EnergyCalculator ecalc, ConfEnergyCalculator confEcalc, EnergyMatrix emat, PruningMatrix pmat, File luteFile) {

		Supplier<ConfSearch> astarFactory = () -> new ConfAStarTree.Builder(emat, pmat)
			.setTraditional()
			.build();

		final File confDBFile = new File("LUTE.conf.db");
		try (ConfDB confdb = new ConfDB(confSpace, confDBFile)) {
			ConfDB.ConfTable confTable = confdb.new ConfTable("lute-test");

			// compute energies for some low-energy confs
			final int numConfs = 100;
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
			//final LUTE.Fitter fitter = LUTE.Fitter.LASSO;
			final LUTE.Fitter fitter = LUTE.Fitter.OLSCG;
			final double maxOverfittingScore = 1.2;
			final double maxRMSE = 0.03;

			confEcalc.resetCounters();

			// do LUTE stuff
			LUTE lute = new LUTE(confSpace);
			ConfSampler sampler = new UniformConfSampler(confSpace, pmat, randomSeed);
			/*ConfSampler sampler = new LowEnergyConfSampler(confSpace, randomSeed, pmat, (rcs) ->
				new ConfAStarTree.Builder(emat, rcs)
					.setTraditional()
					.build()
			);*/
			lute.sampleTuplesAndFit(confEcalc, pmat, confTable, sampler, fitter, maxOverfittingScore, maxRMSE);
			lute.reportConfSpaceSize(pmat);

			// compare conf energies
			double err = 0.0;
			LUTEConfEnergyCalculator luteEcalc = new LUTEConfEnergyCalculator(confSpace, ecalc, new LUTEState(lute.getTrainingSystem()));
			for (ConfSearch.EnergiedConf econf : econfs) {
				double luteEnergy = luteEcalc.calcEnergy(econf.getAssignments());
				log("conf %30s   score %9.4f      energy %9.4f   gap %7.4f      LUTE energy %9.4f   diff %7.4f",
					Arrays.toString(econf.getAssignments()),
					econf.getScore(),
					econf.getEnergy(),
					econf.getEnergy() - econf.getScore(),
					luteEnergy,
					luteEnergy - econf.getEnergy()
				);
				err += (luteEnergy - econf.getEnergy())*(luteEnergy - econf.getEnergy());
			}
			err = Math.sqrt(err/numConfs);
			log("top %d RMSE: %.4f", numConfs, err);

			// save the LUTE reults
			lute.save(luteFile);
		}
	}

	private static void astar(SimpleConfSpace confSpace, EnergyCalculator ecalc, ConfEnergyCalculator confEcalc, EnergyMatrix emat, PruningMatrix pmat, File luteFile) {

		// make a lute ecalc
		LUTEConfEnergyCalculator luteEcalc = new LUTEConfEnergyCalculator(confSpace, ecalc, LUTEIO.read(luteFile));

		ConfAStarTree astar = new ConfAStarTree.Builder(null, pmat)
			.setLUTE(luteEcalc)
			.build();

		log("A* confs:");
		for (int i=0; i<100; i++) {

			ConfSearch.ScoredConf conf = astar.nextConf();
			if (conf == null) {
				log("no more confs");
				break;
			}

			log("%9.4f    conf: %s", conf.getScore(), confSpace.formatConf(conf));
		}
	}

	private static void gmec(SimpleConfSpace confSpace, EnergyCalculator ecalc, ConfEnergyCalculator confEcalc, EnergyMatrix emat, PruningMatrix pmat, File luteFile) {

		// make a lute ecalc
		LUTEConfEnergyCalculator luteEcalc = new LUTEConfEnergyCalculator(confSpace, ecalc, LUTEIO.read(luteFile));

		// do a GMEC search
		Queue.FIFO<ConfSearch.ScoredConf> windowConfs = new LUTEGMECFinder.Builder(pmat, luteEcalc)
			.build()
			.find(1.0);

		log("window confs:");
		while (!windowConfs.isEmpty()) {
			ConfSearch.ScoredConf conf = windowConfs.poll();
			log("%9.4f    conf: %s", conf.getScore(), confSpace.formatConf(conf));
		}
	}

	private static void pfunc(SimpleConfSpace confSpace, EnergyCalculator ecalc, ConfEnergyCalculator confEcalc, EnergyMatrix emat, PruningMatrix pmat, File luteFile) {

		// make a lute ecalc
		LUTEConfEnergyCalculator luteEcalc = new LUTEConfEnergyCalculator(confSpace, ecalc, LUTEIO.read(luteFile));

		final double epsilon = 0.1;

		// make A*
		ConfAStarTree astar = new ConfAStarTree.Builder(null, pmat)
			.setLUTE(luteEcalc)
			.build();

		// calc a pfunc
		LUTEPfunc pfunc = new LUTEPfunc(luteEcalc);
		pfunc.setReportProgress(true);
		pfunc.init(astar, astar.getNumConformations(), epsilon);
		pfunc.compute();

		log("LUTE pfunc bounds: %s", pfunc.makeResult());

		// check the pfunc using traditional methods
		GradientDescentPfunc tpfunc = new GradientDescentPfunc(confEcalc);
		tpfunc.setReportProgress(true);
		tpfunc.init(
			new ConfAStarTree.Builder(emat, pmat)
				.setTraditional()
				.build(),
			new RCs(confSpace).getNumConformations(),
			epsilon
		);
		tpfunc.compute();
		log("traditional pfunc bounds: %s", tpfunc.makeResult());
	}
}
