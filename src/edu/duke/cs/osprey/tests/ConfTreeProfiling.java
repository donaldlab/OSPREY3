package edu.duke.cs.osprey.tests;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.ParallelConfTree;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Stopwatch;

public class ConfTreeProfiling {
	
	public static void main(String[] args)
	throws Exception {
		
		// check the cwd
		String path = new File("").getAbsolutePath();
		if (!path.endsWith("test/DAGK")) {
			throw new Error("This profiler was designed to run in the test/DAGK folder\n\tcwd: " + path);
		}

		// load configuration
		ConfigFileParser cfp = new ConfigFileParser(new String[] {"-c", "KStar.cfg"});
		cfp.loadData();
		
		// multi-thread the energy function
		MultiTermEnergyFunction.setNumThreads(4);
		
		// init a conf space with lots of flexible residues, but no mutations
		final int NumFlexible = 27;
		ArrayList<String> flexRes = new ArrayList<>();
		ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
		for (int i=0; i<NumFlexible; i++) {
			flexRes.add(Integer.toString(i + 1));
			allowedAAs.add(new ArrayList<String>());
		}
		boolean addWt = true;
		boolean doMinimize = false;
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
			flexRes, allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion,
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots
		);
		
		// compute the energy matrix
		System.out.println("\nComputing energy matrix...");
		EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(search.confSpace, search.shellResidues, useERef, addResEntropy);
		emCalc.calcPEM();
		search.emat = emCalc.getEMatrix();
		
		// don't bother with pruning, set all to unpruned
		search.pruneMat = new PruningMatrix(search.confSpace, search.emat.getPruningInterval());
		
		// init the conformation search
		//ConfTree tree = new ConfTree(search);
		ConfTree tree = new ParallelConfTree(search, 2);
		
		// notation below (trialN values in milliseconds):
		// numFlexPos: [trial1, trial2, trial2]
		
		// 2016-05-04
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
		
		System.out.println("\nFinding GMEC among " + tree.getNumConformations().floatValue() + " conformations ...");
		Stopwatch.start();
		tree.nextConf();
		Stopwatch.stop();
		System.out.println("finished in " + Stopwatch.getTime(TimeUnit.MILLISECONDS));
		
		// TODO: check for accuracy
		// energy should be -263.80652436696056
	}
}
