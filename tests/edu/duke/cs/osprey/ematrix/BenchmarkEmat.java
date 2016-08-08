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
		
		String aaNames = "ALA VAL LEU ILE";
		//String aaNames = "ALA";
		//String mutRes = "39";
		String mutRes = "39 43";
		String flexRes = "40 41 42 44 45";
		//String flexRes = "40 41";
		//String flexRes = "40";
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
		
		System.out.println("\nCalculating reference emat...");
		Stopwatch baseStopwatch = new Stopwatch().start();
		EnergyMatrixCalculator emcalc = new EnergyMatrixCalculator(search.confSpace, search.shellResidues, useERef, addResEntropy);
		emcalc.calcPEM();
		search.emat = emcalc.getEMatrix();
		baseStopwatch.stop();
		System.out.println("finished in " + baseStopwatch.getTime());
		
		int[] numThreadsList = { 1, 2, 4, 8 };
		
		// what energy function generator to use?
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		
		/* NOTE: the gpu is pretty slow at the very small energy functions (eg pairwise residues) compared to the cpu
		int maxNumThreads = numThreadsList[numThreadsList.length - 1];
		GpuQueuePool pool = new GpuQueuePool(maxNumThreads, 1);
		GpuQueuePool pool = new GpuQueuePool(1, maxNumThreads);
		EnergyFunctionGenerator egen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), pool);
		*/
		
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, search.confSpace, search.shellResidues);
		
		System.out.println("\nBenchmarking Emat calculation, main thread...");
		Stopwatch mainStopwatch = new Stopwatch().start();
		EnergyMatrix emat = new SimpleEnergyMatrixCalculator(ecalc).calcEnergyMatrix();
		mainStopwatch.stop();
		System.out.println(String.format("Speedup: %.2fx", (float)baseStopwatch.getTimeNs()/mainStopwatch.getTimeNs()));
		checkEmat(search.emat, emat);
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		
		for (int numThreads : numThreadsList) {
			
			tasks.start(numThreads);
			
			System.out.println("\nBenchmarking Emat calculation, " + numThreads + " task thread(s)...");
			Stopwatch taskStopwatch = new Stopwatch().start();
			emat = new SimpleEnergyMatrixCalculator(ecalc).calcEnergyMatrix(tasks);
			taskStopwatch.stop();
			System.out.println(String.format("Speedup: %.2fx", (float)baseStopwatch.getTimeNs()/taskStopwatch.getTimeNs()));
			checkEmat(search.emat, emat);
		
			tasks.stopAndWait(10000);
		}
		
		//pool.cleanup();
	}

	private static void checkEmat(EnergyMatrix exp, EnergyMatrix obs) {
		
		for (int pos1=0; pos1<exp.getNumPos(); pos1++) {
			for (int rc1=0; rc1<exp.getNumConfAtPos(pos1); rc1++) {
			
				// singles
				checkEnergy(exp.getOneBody(pos1, rc1), obs.getOneBody(pos1, rc1));
				
				// pairs
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<exp.getNumConfAtPos(pos2); rc2++) {
						
						checkEnergy(exp.getPairwise(pos1, rc1, pos2, rc2), obs.getPairwise(pos1, rc1, pos2, rc2));
					}
				}
			}
		}
	}

	private static void checkEnergy(double exp, double obs) {
		
		// minimized energies aren't super precise
		final double Epsilon = 1e-3;
		
		double absErr = Math.abs(exp - obs);
		double relErr = absErr/Math.abs(exp);
		if (relErr > Epsilon) {
			System.out.println(String.format("\tWARNING: low energy precision!  expected: %12.6f  observed: %12.6f  absErr: %12.6f  relErr: %12.6f",
				exp, obs, absErr, relErr
			));
		}
	}
}
