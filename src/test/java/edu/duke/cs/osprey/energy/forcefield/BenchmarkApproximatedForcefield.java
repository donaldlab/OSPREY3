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

package edu.duke.cs.osprey.energy.forcefield;


import edu.duke.cs.osprey.Benchmark;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.approximation.ApproximatorMatrix;
import edu.duke.cs.osprey.energy.approximation.ApproximatorMatrixCalculator;
import edu.duke.cs.osprey.energy.approximation.ResidueInteractionsApproximator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Progress;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.util.function.BiConsumer;
import java.util.function.BiFunction;

import static edu.duke.cs.osprey.tools.Log.log;


public class BenchmarkApproximatedForcefield {

	public static void main(String[] args) {

		// make a conf space
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A68").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // arg
		strand.flexibility.get("A6" ).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // hie
		strand.flexibility.get("A24").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // lys
		strand.flexibility.get("A32").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // asp
		strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // glu
		strand.flexibility.get("A15").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // cys
		strand.flexibility.get("A17").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // gly
		strand.flexibility.get("A31").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // pro
		strand.flexibility.get("A2" ).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // ala
		strand.flexibility.get("A11").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // val
		strand.flexibility.get("A36").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // ile
		strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // leu
		strand.flexibility.get("A13").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // met
		strand.flexibility.get("A55").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // phe
		strand.flexibility.get("A48").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // tyr
		// no trp
		strand.flexibility.get("A69").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // ser
		strand.flexibility.get("A14").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // thr
		strand.flexibility.get("A23").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // asn
		strand.flexibility.get("A43").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // gln

		/*
		for (String resNum : Arrays.asList("A16", "A17", "A18", "A19", "A20", "A21")) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}
		*/
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		BiFunction<ConfEnergyCalculator,String,EnergyMatrix> calcEmat = (confEcalc, label) ->
			new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				//.setTripleCorrectionThreshold(10.0)
				.setCacheFile(new File("benchmark." + label + ".emat"))
				.build()
				.calcEnergyMatrix();

		BiConsumer<EnergyMatrix,ConfEnergyCalculator> calcConfs = (emat, confEcalc) -> {
			ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
				.setTraditional()
				.build();
			Progress progress = new Progress(100);
			for (int i=0; i<progress.getTotalWork(); i++) {
				confEcalc.calcEnergy(astar.nextConf());
				progress.incrementProgress();
			}
		};

		// get an energy calculator
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build()
		) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				//.setEnergyPartition(EnergyPartition.Traditional)
				.setEnergyPartition(EnergyPartition.AllOnPairs)
				.build();

			// calc the approximator matrix
			ApproximatorMatrix amat = new ApproximatorMatrixCalculator(confEcalc)
				.setCacheFile(new File("benchmark.amat"))
				.calc();

			ConfEnergyCalculator confEcalcApprox = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				.setEnergyPartition(confEcalc.epart)
				.setApproximatorMatrix(amat)
				.setApproximationErrorBudget(1e-1)
				//.setApproximationErrorBudget(10)
				//.setApproximationErrorBudget(Double.POSITIVE_INFINITY)
				.build();

			// calc energy matrices

			/* TEMP
			log("\ncalculating real emat...");
			Stopwatch realEmatStopwatch = new Stopwatch().start();
			EnergyMatrix emat = calcEmat.apply(confEcalc, "real");
			log("\tdone in %s", realEmatStopwatch.stop().getTime(2));

			log("\ncalculating approx emat...");
			Stopwatch approxEmatStopwatch = new Stopwatch().start();
			EnergyMatrix ematApprox = calcEmat.apply(confEcalcApprox, "approx");
			log("\tdone in %s", approxEmatStopwatch.stop().getTime(2));
			*/

			// calc full conf minimizations

			/*
			log("\ncalculating real confs...");
			Stopwatch realConfsStopwatch = new Stopwatch().start();
			calcConfs.accept(emat, confEcalc);
			log("\tdone in %s", realConfsStopwatch.stop().getTime(2));

			log("\ncalculating approx confs...");
			Stopwatch approxConfsStopwatch = new Stopwatch().start();
			calcConfs.accept(ematApprox, confEcalcApprox);
			log("\tdone in %s", approxConfsStopwatch.stop().getTime(2));
			*/

			// benchmark minimizing a full conf (with the wild-type conf)
			RCTuple tuple = new RCTuple();
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				tuple.pos.add(pos.index);
				tuple.RCs.add(pos.resConfs.size() - 1);
			}

			{
				ResidueInteractions inters = EnergyPartition.makeFragment(confSpace, null, false, tuple);
				double energyReal = confEcalc.calcEnergy(tuple, inters).energy;
				double energyApprox = confEcalcApprox.calcEnergy(tuple, inters).energy;

				ResidueInteractionsApproximator approximator = amat.get(tuple, inters, confEcalcApprox.approximationErrorBudget);

				Benchmark bmMol = new Benchmark(1, 100, 1000, () -> confSpace.makeDiscreteMolecule(tuple));
				Benchmark bmReal = new Benchmark(1, 10, 20, () -> confEcalc.calcEnergy(tuple, inters));
				Benchmark bmApprox = new Benchmark(1, 10, 20, () -> confEcalcApprox.calcEnergy(tuple, inters));

				log("mol=[%8.2f ms   %6.1f ops]   real=[%8.4f   %8.2f ms   %6.1f ops]   approx=[%8.4f   %8.2f ms   %6.1f ops]   error=%8.4f   inters=%3d/%3d (%5.1f%%)   speedup=%4.1fx",
					bmMol.stopwatch.getTimeMs(), bmMol.opsPerSecond,
					energyReal, bmReal.stopwatch.getTimeMs(), bmReal.opsPerSecond,
					energyApprox, bmApprox.stopwatch.getTimeMs(), bmApprox.opsPerSecond,
					Math.abs(energyReal - energyApprox),
					approximator.approxInters.size(), inters.size(), 100f*approximator.approxInters.size()/inters.size(),
					bmApprox.opsPerSecond/bmReal.opsPerSecond
				);
			}

			/*
			// breakdown each residue
			for (SimpleConfSpace.Position pos : confSpace.positions) {
			//{ SimpleConfSpace.Position pos = confSpace.positions.get(6);

				// get the wild type residue
				SimpleConfSpace.ResidueConf rc = pos.resConfs.get(pos.resConfs.size() - 1);
				assert (rc.type == SimpleConfSpace.ResidueConf.Type.WildType);

				RCTuple tuple = new RCTuple(pos.index, rc.index);

				// only shell interactions for now
				ResidueInteractions inters = ResInterGen.of(confSpace)
					.addIntra(pos.index)
					.addShell(pos.index)
					.make();

				double energyReal = confEcalc.calcEnergy(tuple, inters).energy;
				double energyApprox = confEcalcApprox.calcEnergy(tuple, inters).energy;

				ResidueInteractionsApproximator approximator = amat.get(tuple, inters, confEcalcApprox.approximationErrorBudget);

				//int numRuns = 1000;
				int numRuns = 1000;
				Benchmark bmMol = new Benchmark(1, 100, numRuns, () -> confSpace.makeDiscreteMolecule(tuple));
				Benchmark bmReal = new Benchmark(1, 100, numRuns, () -> confEcalc.calcEnergy(tuple, inters));
				Benchmark bmApprox = new Benchmark(1, 100, numRuns, () -> confEcalcApprox.calcEnergy(tuple, inters));

				log("%3s %3s %2d:%2d   mol=[%8.2f ms   %6.1f ops]   real=[%8.4f   %8.2f ms   %6.1f ops]   approx=[%8.4f   %8.2f ms   %6.1f ops]   error=%8.4f   inters=%3d/%3d (%5.1f%%)   speedup=%4.1fx",
					pos.resNum, pos.resTypes.get(0),
					pos.index, rc.index,
					bmMol.stopwatch.getTimeMs(), bmMol.opsPerSecond,
					energyReal, bmReal.stopwatch.getTimeMs(), bmReal.opsPerSecond,
					energyApprox, bmApprox.stopwatch.getTimeMs(), bmApprox.opsPerSecond,
					Math.abs(energyReal - energyApprox),
					approximator.inters.size(), inters.size(), 100f*approximator.inters.size()/inters.size(),
					bmApprox.opsPerSecond/bmReal.opsPerSecond
				);
			}
			*/
		}
	}
}
