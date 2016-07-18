package edu.duke.cs.osprey.gpu;

import java.util.ArrayList;
import java.util.Arrays;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.tools.Stopwatch;

public class BenchmarkForcefieldKernel extends TestBase {
	
	public static void main(String[] args)
	throws Exception {
		
		initDefaultEnvironment();
		
		// make a search problem
		System.out.println("Building search problem...");
		
		String aaNames = "ALA VAL LEU ILE";
		//String aaNames = "ALA";
		String mutRes = "39";
		//String mutRes = "39 43";
		//String flexRes = "";
		String flexRes = "40";
		//String flexRes = "40 41";
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
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots
		);
		
		benchmarkEfunc(search);
		//benchmarkEmat(search);
	}
	
	private static void benchmarkEmat(SearchProblem search)
	throws Exception {
		// TODO
	}
	
	private static void benchmarkEfunc(SearchProblem search)
	throws Exception {
		
		int numRuns = 1000;
		double energy = 0;
		
		// build the energy functions
		System.out.println("Building energy functions...");
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		ForcefieldParams ffparams = makeDefaultFFParams();
		GpuEnergyFunctionGenerator gpuegen = new GpuEnergyFunctionGenerator(ffparams);
		EnergyFunction efunc = egen.fullConfEnergy(search.confSpace, search.shellResidues);
		GpuForcefieldEnergy gpuefunc = gpuegen.fullConfEnergy(search.confSpace, search.shellResidues);
		System.out.println("\tnum atom pairs: " + gpuefunc.getForcefieldEnergy().getNumAtomPairs());
		System.out.println(String.format("\tGPU memory used: %.2f Mib", (double)gpuefunc.getGpuBytesNeeded()/1024/1024));
		
		// benchmark the cpu
		System.out.println("\nBenchmarking CPU...");
		Stopwatch cpuStopwatch = new Stopwatch().start(); 
		for (int i=0; i<numRuns; i++) {
			energy = efunc.getEnergy();
		}
		cpuStopwatch.stop();
		System.out.println("\tfinished in " + cpuStopwatch.getTime(2));
		System.out.println("\tenergy: " + energy);
		
		// benchmark the gpu
		System.out.println("\nBenchmarking GPU...");
		Stopwatch gpuStopwatch = new Stopwatch().start();
		for (int i=0; i<numRuns; i++) {
			energy = gpuefunc.getEnergy();
		}
		gpuStopwatch.stop();
		gpuefunc.cleanup();
		System.out.println("\tfinished in " + gpuStopwatch.getTime(2));
		System.out.println("\tenergy: " + energy);
		
		// calculate speedup
		System.out.println(String.format("\nspeedup: %.2fx", (double)cpuStopwatch.getTimeNs()/gpuStopwatch.getTimeNs())); 
	}
}
