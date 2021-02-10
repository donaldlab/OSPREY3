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

package edu.duke.cs.osprey;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.TimeUnit;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class ConfTreeProfiling {
	
	public static void main(String[] args)
	throws Exception {
		
		// check the cwd
		String path = new File("").getAbsolutePath();
		if (!path.endsWith("examples/DAGK")) {
			throw new Error("This profiler was designed to run in the examples/DAGK folder\n\tcwd: " + path);
		}

		// load configuration
		ConfigFileParser cfp = new ConfigFileParser();
		cfp.loadData();
		
		// multi-thread the energy function
		MultiTermEnergyFunction.setNumThreads(4);
		
		// init a conf space with lots of flexible residues, but no mutations
		
		boolean doMinimize = true;
		//final int NumFlexible = 16;
		//final int NumFlexible = 25;
		final int NumFlexible = 40;
		
		//boolean doMinimize = false;
		//final int NumFlexible = 27;
		//final int NumFlexible = 34;
		//final int NumFlexible = 55;
		
		ArrayList<String> flexRes = new ArrayList<>();
		ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
		for (int i=0; i<NumFlexible; i++) {
			flexRes.add(Integer.toString(i + 1));
			allowedAAs.add(new ArrayList<String>());
		}
		boolean addWt = true;
		boolean useEpic = false;
		boolean useTupleExpansion = false;
		boolean useEllipses = false;
		boolean useERef = false;
		boolean addResEntropy = false;
		boolean addWtRots = true;
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		SearchProblem search = new SearchProblem(
			"energyMatrixProfiling",
			"2KDC.P.forOsprey.pdb", 
			flexRes, allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null, 
                        false, new ArrayList<>()
		);
		
		// compute the energy matrix
		File ematFile = new File(String.format("emat.%s%d.dat", doMinimize ? "min." : "", NumFlexible));
		if (ematFile.exists()) {
			System.out.println("\nReading energy matrix...");
			search.emat = (EnergyMatrix)ObjectIO.readObject(ematFile.getAbsolutePath(), true);
		}
		if (search.emat == null) {
			System.out.println("\nComputing energy matrix...");
			EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(search.confSpace, search.shellResidues, useERef, addResEntropy);
			emCalc.calcPEM();
			search.emat = emCalc.getEMatrix();
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
		}
		
		// don't bother with pruning, set all to unpruned
		search.pruneMat = new PruningMatrix(search.confSpace, search.emat.getPruningInterval());
		
		/* TEMP: set some positions to just one residue, to test A* order
		rcs.set(0, new int[] { 0 });
		rcs.set(1, new int[] { 0 });
		rcs.set(2, new int[] { 7 });
		rcs.set(3, new int[] { 16 });
		rcs.set(4, new int[] { 16 });
		rcs.set(5, new int[] { 0 });
		*/
		
		// init the conformation search
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setMPLP(new ConfAStarTree.MPLPBuilder()
				.setNumIterations(5)
			).setShowProgress(true)
			.build();
		
		// notation below (trialN values in milliseconds):
		// numFlexPos: [trial1, trial2, trial2]
		
		// 2016-05-03
		// BEFORE OPTIMIZATIONS
		// 27: [36503, 37969, 36664]
		
		// after roughly 2x energy matrix read speedup
		// 27: [25663, 25565, 25646] => 1.45x speedup over benchmark
		
		// optimize ConfTree a bit
		// 27: [18446, 18387, 18470] => 2.01x speedup over benchmark
		
		// implement lazy instantiation of higher-order terms
		// 27: [12963, 13017, 12885] => 2.86x speedup over benchmark
		
		// implement parallel expansion of nodes (2 threads)
		// 27x2: [14847, 15117, 15544] => rather disappointing, really
		// this workload must be memory-bound, so parallel computations hurt more than help
		
		// NB: turning off dynamic A* bumped the run time to at least 5 minutes
		// I stopped waiting after that
		// dynamic A* makes a HUGE difference!!
		
		// 2016-05-04
		// another day, another something... dunno what it is. re-doing benchmarks
		
		// after 2x energy matrix speedup
		// 27: [31356, 31728, 31215]
		
		// Bah. I'm lazy and don't want to re-benchmark the before-optimization state
		// let's just assume this weird daily thing is linear and say we're 1.23x slower today
		// so: pre-opt should be:
		// 27: [44777, 46575, 44975]
		// which means the energy matrix improvements are still about a 1.45x speedup
		// so that checks out
		
		// yesterday's ConfTree optimizations:
		// 27: [15723, 15778, 15897] => 2.88x speedup over benchmark
		
		// today's ConfTree optimizations:
		// 27: [2619, 2665, 2663] => 17.15x speedup over benchmark!! =D
		
		// NOTE: addWTRots bug fix on 2016-05-05 changed energy matrix values!
		// so all benchmarks after that are incomparable to benchmarks before it
		// also, the newer progress reporting does impose a small performance penalty
		
		// sooo..... let's start some new benchmarks
		// 2016-05-13
		// current state of code:
		// 34:   [22759, 22687, 22572]
		
		// after minor optimizations, haven't dropped the big guns just yet... =P
		// 34:   [18846, 18921, 18962] => 1.20x speedup over benchmark
		
		// 2016-05-14 (didn't bother re-benchmarking today)
		// after differential score calculations
		// 34:   [1325, 1326, 1337] => 17.06x speedup over benchmark!! =D
		
		// this test run is too short now... need something longer
		// 55:   [24873, 24501, 25076]
		
		// 2016-05-15 (didn't bother re-benchmarking again)
		// after more minor optimizations, fixed memory usage
		// 55:   [19785, 20181, 20118] => 1.24x speedup over benchmark
		
		// optimize for memory usage at the nodes, hopefully this won't slow down too much
		// 55:   [20082, 20240, 20227] => about the same, not bad at all! =)
		
		System.out.println("\nFinding GMEC among " + tree.getNumConformations().doubleValue() + " conformations ...");
		Stopwatch stopwatch = new Stopwatch();
		stopwatch.start();
		ConfSearch.ScoredConf conf = tree.nextConf();
		stopwatch.stop();
		System.out.println("finished in " + stopwatch.getTime(TimeUnit.MILLISECONDS));
		System.out.println("conf:     " + Arrays.toString(conf.getAssignments()));
		
		double observedEnergy = conf.getScore();
		System.out.println(String.format("energy: %16.10f", observedEnergy));
		
		if (doMinimize) {
			
			// make sure we still have the right answer!
			Map<Integer,int[]> expectedConfs = new TreeMap<>();
			expectedConfs.put(16, new int[] { 0, 0, 7, 16, 16, 0, 1, 0, 25, 6, 0, 21, 0, 0, 0, 6 });
			expectedConfs.put(25, new int[] { 0, 0, 7, 16, 16, 0, 1, 0, 25, 6, 0, 21, 0, 0, 0, 6, 3, 2, 10, 1, 0, 0, 0, 1, 0 });
			expectedConfs.put(40, new int[] { 0, 0, 7, 16, 16, 0, 1, 0, 25, 6, 0, 21, 0, 0, 0, 6, 3, 2, 10, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 3 });
			checkConf(expectedConfs.get(NumFlexible), conf.getAssignments());
			
			// make sure the energy matches
			Map<Integer,Double> expectedEnergies = new TreeMap<>();
			expectedEnergies.put(16, -160.9788154538);
			expectedEnergies.put(25, -261.3928439504);
			expectedEnergies.put(40, -421.9530270377);
			checkEnergy(expectedEnergies.get(NumFlexible), observedEnergy);
			
		} else {
		
			// make sure we still have the right answer!
			Map<Integer,int[]> expectedConfs = new TreeMap<>();
			expectedConfs.put(27, new int[] { 0, 6, 7, 0, 16, 0, 0, 6, 25, 6, 0, 0, 0, 0, 0, 0, 16, 2, 12, 1, 0, 15, 0, 1, 0, 0, 0 });
			expectedConfs.put(34, new int[] { 0, 6, 7, 0, 16, 0, 0, 6, 25, 6, 0, 0, 0, 0, 0, 0, 16, 2, 12, 1, 0, 15, 0, 1, 0, 0, 0, 0, 0, 0, 0, 29, 0, 0 });
			expectedConfs.put(55, new int[] { 0, 6, 7, 0, 16, 0, 0, 6, 25, 6, 0, 0, 0, 0, 0, 0, 16, 2, 12, 1, 0, 15, 0, 1, 0, 0, 0, 0, 0, 0, 0, 29, 0, 0, 1, 2, 1, 0, 0, 3, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
			checkConf(expectedConfs.get(NumFlexible), conf.getAssignments());
			
			// make sure the energy matches
			Map<Integer,Double> expectedEnergies = new TreeMap<>();
			expectedEnergies.put(27, -260.91555715297517);
			expectedEnergies.put(34, -346.32024675046176);
			expectedEnergies.put(55, -514.1055956242977);
			checkEnergy(expectedEnergies.get(NumFlexible), observedEnergy);
		}
	}
	
	private static void checkConf(int[] expected, int[] observed) {
		if (!Arrays.equals(expected, observed)) {
			System.out.println("expected: " + Arrays.toString(expected));
			throw new Error("GMEC changed! Undo that \"optimization!\"");
		}
	}
	
	private static void checkEnergy(double expected, double observed) {
		final double Epsilon = 1e-14;
		double relErr = Math.abs(expected - observed)/expected;
		if (relErr > Epsilon) {
			throw new Error(String.format(
				"Energy is wrong, undo that 'optimization'!\n\texpected: %.15f\n\tcalculated: %.15f\n\trelErr: %.15f",
				expected, observed, relErr
			));
		} else {
			System.out.println(String.format("Energy is correct\n\trelErr: %.15f", relErr));
		}
	}
}
