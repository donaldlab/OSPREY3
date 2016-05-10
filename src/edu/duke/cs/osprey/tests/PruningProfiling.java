package edu.duke.cs.osprey.tests;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;

public class PruningProfiling {
	
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
		// n=300 gives about 2.5e241 conformations to prune
		final int NumFlexible = 300;
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
		
		// compute/read the energy matrix
		File ematFile = new File(String.format("emat.%d.dat", NumFlexible));
		if (ematFile.exists()) {
			System.out.println("\nReading energy matrix...");
			search.emat = (EnergyMatrix)ObjectIO.readObject(ematFile.getAbsolutePath(), true);
		} else {
			System.out.println("\nComputing energy matrix...");
			EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(search.confSpace, search.shellResidues, useERef, addResEntropy);
			emCalc.calcPEM();
			search.emat = emCalc.getEMatrix();
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
		}
		
		// init pruning
		double pruningInterval = 0;
		search.pruneMat = new PruningMatrix(search.confSpace, pruningInterval);
		boolean typeDep = false;
		double boundsThreshold = 100; // config default
		int algoOption = 3; // these two kind of do the same thing: turn on pair pruning
		boolean useFlags = true; //
		boolean useTriples = false;
		boolean preDacs = false;
		double stericThreshold = 100; // config default
		PruningControl pruner = new PruningControl(
			search, pruningInterval, typeDep, boundsThreshold, algoOption, 
			useFlags, useTriples, preDacs, useEpic, useTupleExpansion, stericThreshold);
		
		// notation below (trialN values in milliseconds):
		// numFlexPos: [trial1, trial2, trial2]
		
		// 2016-05-10
		// BEFORE OPTIMIZATIONS
		// 300: [18510, 18428, 18635]
		
		System.out.println("\nPruning " + search.confSpace.getNumConformations().doubleValue() + " conformations...");
		Stopwatch.start();
		pruner.prune();            
		Stopwatch.stop();
		System.out.println("finished in " + Stopwatch.getTime(TimeUnit.MILLISECONDS));
	}
}
