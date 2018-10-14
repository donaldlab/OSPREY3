/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

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
import edu.duke.cs.osprey.kstar.BBKStar;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;


public class LUTELab {

	public static void main(String[] args) {
		//gmecMain();
		kstarMain();
	}

	public static void kstarMain() {

		//TestKStar.ConfSpaces confSpaces = TestKStar.make1GUA11();
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();

		//train("protein", confSpaces.protein);
		//train("ligand", confSpaces.ligand);
		//train("complex", confSpaces.complex);
		//kstarPfunc(confSpaces);
		//kstar(confSpaces);
		bbkstar(confSpaces);
	}

	private static void train(String name, SimpleConfSpace confSpace) {

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			//.setParallelism(Parallelism.makeCpu(8))
			.setParallelism(Parallelism.make(8, 3, 32))
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
					.build()
					.calcReferenceEnergies()
				).build();

			// compute energy matrix
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(new File(String.format("LUTE.%s.emat.dat", name)))
				.build()
				.calcEnergyMatrix();

			// run DEE (super important for good LUTE fits!!)
			PruningMatrix pmat = new SimpleDEE.Runner()
				.setSinglesThreshold(100.0)
				.setPairsThreshold(100.0)
				.setGoldsteinDiffThreshold(10.0)
				.setShowProgress(true)
				.setCacheFile(new File(String.format("LUTE.%s.pmat.dat", name)))
				.setParallelism(Parallelism.makeCpu(8))
				.run(confSpace, emat);

			final File confDBFile = new File(String.format("LUTE.%s.conf.db", name));
			try (ConfDB confdb = new ConfDB(confSpace, confDBFile)) {
				ConfDB.ConfTable confTable = confdb.new ConfTable("lute");

				log("\nLUTE:\n");

				final int randomSeed = 12345;
				final LUTE.Fitter fitter = LUTE.Fitter.OLSCG;
				final double maxOverfittingScore = 1.5;
				final double maxRMSE = 0.1;

				confEcalc.resetCounters();

				// compute LUTE fit
				LUTE lute = new LUTE(confSpace);
				ConfSampler sampler = new RandomizedDFSConfSampler(confSpace, pmat, randomSeed);
				lute.sampleTuplesAndFit(confEcalc, emat, pmat, confTable, sampler, fitter, maxOverfittingScore, maxRMSE);
				lute.reportConfSpaceSize(pmat);
				lute.save(new File(String.format("LUTE.%s.dat", name)));
			}
		}
	}

	private static void kstarPfunc(TestKStar.ConfSpaces confSpaces) {

		// calculate the pfunc of the complex sequences
		SimpleConfSpace confSpace = confSpaces.complex;
		String name = "complex";

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces.complex, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(8))
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(new File(String.format("LUTE.%s.emat.dat", name)))
				.build()
				.calcEnergyMatrix();

			PruningMatrix pmat = new SimpleDEE.Runner()
				.setSinglesThreshold(100.0)
				.setPairsThreshold(100.0)
				.setGoldsteinDiffThreshold(10.0)
				.setShowProgress(true)
				.setCacheFile(new File(String.format("LUTE.%s.pmat.dat", name)))
				.setParallelism(Parallelism.makeCpu(8))
				.run(confSpace, emat);

			LUTEState luteState = LUTEIO.read(new File(String.format("LUTE.%s.dat", name)));
			LUTEConfEnergyCalculator luteEcalc = new LUTEConfEnergyCalculator(confSpaces.complex, luteState);

			// enumerate all the sequences
			List<Sequence> sequences = new ArrayList<>();
			for (List<SimpleConfSpace.Position> mutablePositions : MathTools.powersetUpTo(confSpace.positions, 1)) {

				// collect the mutations (res types except for wild type) for these positions into a simple list list
				List<List<String>> resTypes = new ArrayList<>();
				for (SimpleConfSpace.Position pos : mutablePositions) {
					resTypes.add(pos.resFlex.resTypes.stream()
						.filter((resType) -> !resType.equals(pos.resFlex.wildType))
						.collect(Collectors.toList())
					);
				}

				// enumerate all the combinations of res types
				for (List<String> mutations : MathTools.cartesianProduct(resTypes)) {
					Sequence sequence = confSpace.makeWildTypeSequence();
					for (int i=0; i<mutablePositions.size(); i++) {
						sequence.set(mutablePositions.get(i).resNum, mutations.get(i));
					}
					sequences.add(sequence);
				}
			}

			for (Sequence sequence : sequences) {

				RCs unprunedRCs = sequence.makeRCs(confSpace);
				RCs prunedRCs = new RCs(unprunedRCs, pmat);

				ConfAStarTree astar = new ConfAStarTree.Builder(null, prunedRCs)
					.setLUTE(luteEcalc)
					.build();

				LUTEPfunc pfunc = new LUTEPfunc(luteEcalc);
				pfunc.setReportProgress(false);
				pfunc.init(astar, unprunedRCs.getNumConformations(), 0.01);
				pfunc.compute();
				log("LUTE pfunc bounds for %s: %s", sequence.toString(Sequence.Renderer.ResTypeMutations), pfunc.makeResult());
			}
		}
	}

	private static void kstar(TestKStar.ConfSpaces confSpaces) {

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces.complex, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(8))
			.build()) {

			KStar.Settings settings = new KStar.Settings.Builder()
				.setEpsilon(0.01)
				.setStabilityThreshold(null)
				.setMaxSimultaneousMutations(1)
				.addScoreConsoleWriter()
				.build();
			KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, settings);
			for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

				// how should we define energies of conformations?
				LUTEState luteState = LUTEIO.read(new File(String.format("LUTE.%s.dat", info.id)));
				LUTEConfEnergyCalculator luteConfEcalc = new LUTEConfEnergyCalculator(info.confSpace, luteState);
				info.confEcalc = luteConfEcalc;

				// load the pruning matrix
				PruningMatrix pmat = SimpleDEE.read(info.confSpace, new File(String.format("LUTE.%s.pmat.dat", info.id)));

				// how should confs be ordered and searched?
				info.confSearchFactory = (rcs) -> {
					rcs = new RCs(rcs, pmat);
					return new ConfAStarTree.Builder(null, rcs)
						.setLUTE(luteConfEcalc)
						.build();
				};
			}
			kstar.run();
		}
	}

	private static void bbkstar(TestKStar.ConfSpaces confSpaces) {

		try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(8))
			.build()) {

			KStar.Settings kstarSettings = new KStar.Settings.Builder()
				.setEpsilon(0.01)
				.setStabilityThreshold(null)
				.setMaxSimultaneousMutations(1)
				.addScoreConsoleWriter()
				.build();
			BBKStar.Settings bbkstarSettings = new BBKStar.Settings.Builder()
				.setNumBestSequences(25)
				.build();
			BBKStar bbkstar = new BBKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, kstarSettings, bbkstarSettings);
			for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {

				// how should we define energies of conformations?
				LUTEState luteState = LUTEIO.read(new File(String.format("LUTE.%s.dat", info.id)));
				LUTEConfEnergyCalculator luteConfEcalc = new LUTEConfEnergyCalculator(info.confSpace, luteState);
				info.confEcalcMinimized = luteConfEcalc;

				// load the pruning matrix
				PruningMatrix pmat = SimpleDEE.read(info.confSpace, new File(String.format("LUTE.%s.pmat.dat", info.id)));

				// how should confs be ordered and searched?
				info.confSearchFactoryMinimized = (rcs) -> {
					rcs = new RCs(rcs, pmat);
					return new ConfAStarTree.Builder(null, rcs)
						.setLUTE(luteConfEcalc)
						.build();
				};

				// for rigid energies, just use the minimized energies, since they're all the same for LUTE
				info.confSearchFactoryRigid = info.confSearchFactoryMinimized;
			}
			bbkstar.run();
		}
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
			ConfSampler sampler = new RandomizedDFSConfSampler(confSpace, pmat, randomSeed);
			/*ConfSampler sampler = new LowEnergyConfSampler(confSpace, randomSeed, pmat, (rcs) ->
				new ConfAStarTree.Builder(emat, rcs)
					.setTraditional()
					.build()
			);*/
			lute.sampleTuplesAndFit(confEcalc, emat, pmat, confTable, sampler, fitter, maxOverfittingScore, maxRMSE);
			lute.reportConfSpaceSize(pmat);

			// compare conf energies
			double err = 0.0;
			LUTEConfEnergyCalculator luteEcalc = new LUTEConfEnergyCalculator(confSpace, new LUTEState(lute.getTrainingSystem()));
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
		LUTEConfEnergyCalculator luteEcalc = new LUTEConfEnergyCalculator(confSpace, LUTEIO.read(luteFile));

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
		LUTEConfEnergyCalculator luteEcalc = new LUTEConfEnergyCalculator(confSpace, LUTEIO.read(luteFile));

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
		LUTEConfEnergyCalculator luteEcalc = new LUTEConfEnergyCalculator(confSpace, LUTEIO.read(luteFile));

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
