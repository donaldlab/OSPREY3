package edu.duke.cs.osprey.ematrix;

import java.util.ArrayList;
import java.util.Arrays;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.Stopwatch;

public class BenchmarkEmat extends TestBase {
	
	public static void main(String[] args)
	throws Exception {
		
		initDefaultEnvironment();
		
		// don't use energy function-level parallelism
		MultiTermEnergyFunction.setNumThreads(1);
		
		// make a search problem
		System.out.println("Building search problem...");
		
		//String aaNames = "ALA VAL LEU ILE";
		String aaNames = "ALA";
		String mutRes = "39 43";
		//String flexRes = "40 41 42 44 45";
		String flexRes = "40 41";
		ArrayList<String> flexResList = new ArrayList<>();
		ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
		for (String res : mutRes.split(" ")) {
			if (!res.isEmpty()) {
				flexResList.add(res);
				allowedAAs.add(new ArrayList<>(Arrays.asList(aaNames.split(" "))));
			}
		}
		for (String res : flexRes.split(" ")) {
			if (!res.isEmpty()) {
				flexResList.add(res);
				allowedAAs.add(new ArrayList<>());
			}
		}
		boolean doMinimize = true;
		boolean addWt = true;
		boolean useEpic = false;
		boolean useTupleExpansion = false;
		boolean useEllipses = false;
		boolean useERef = false;
		boolean addResEntropy = false;
		boolean addWtRots = false;
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		SearchProblem search = new SearchProblem(
			"test", "test/1CC8/1CC8.ss.pdb", 
			flexResList, allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion,
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null
		);
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		//GpuEnergyFunctionGenerator gpuegen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), new GpuQueuePool(2, 2));
		
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, search.confSpace, search.shellResidues);
		
		System.out.println("\nBenchmarking Emat calculation, 1 thread...");
		Stopwatch thread1Stopwatch = new Stopwatch().start();
		new SimpleEnergyMatrixCalculator(ecalc).calcEnergyMatrix();
		thread1Stopwatch.stop();
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(2);
		
		System.out.println("\nBenchmarking Emat calculation, 2 threads...");
		Stopwatch thread2Stopwatch = new Stopwatch().start();
		new SimpleEnergyMatrixCalculator(ecalc).calcEnergyMatrix(tasks);
		thread2Stopwatch.stop();
		
		tasks.stopAndWait(10000);
	}
}
