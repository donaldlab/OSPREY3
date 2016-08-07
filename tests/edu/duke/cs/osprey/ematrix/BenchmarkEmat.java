package edu.duke.cs.osprey.ematrix;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.gpu.GpuQueuePool;
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
		boolean doMinimize = false;
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
		EnergyMatrixCalculator emcalc = new EnergyMatrixCalculator(search.confSpace, search.shellResidues, useERef, addResEntropy);
		emcalc.calcPEM();
		search.emat = emcalc.getEMatrix();
		
		//EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		EnergyFunctionGenerator egen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), new GpuQueuePool(1, 1));
		// NOTE: the gpu is pretty slow at the very small energy functions (eg pairwise residies) compared to the cpu
		
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, search.confSpace, search.shellResidues);
		
		System.out.println("\nBenchmarking Emat calculation, main thread...");
		Stopwatch thread1Stopwatch = new Stopwatch().start();
		EnergyMatrix emat = new SimpleEnergyMatrixCalculator(ecalc).calcEnergyMatrix();
		checkEmat(search.emat, emat);
		thread1Stopwatch.stop();
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(1);
		
		System.out.println("\nBenchmarking Emat calculation, 1 task thread...");
		Stopwatch tasks1Stopwatch = new Stopwatch().start();
		EnergyMatrix emat1 = new SimpleEnergyMatrixCalculator(ecalc).calcEnergyMatrix(tasks);
		checkEmat(search.emat, emat1);
		tasks1Stopwatch.stop();
		System.out.println(String.format("Speedup: %.2fx", (float)thread1Stopwatch.getTimeNs()/tasks1Stopwatch.getTimeNs()));
		
		tasks.stopAndWait(10000);
		tasks.start(2);
		
		System.out.println("\nBenchmarking Emat calculation, 2 task threads...");
		Stopwatch tasks2Stopwatch = new Stopwatch().start();
		EnergyMatrix emat2 = new SimpleEnergyMatrixCalculator(ecalc).calcEnergyMatrix(tasks);
		checkEmat(search.emat, emat2);
		tasks2Stopwatch.stop();
		System.out.println(String.format("Speedup: %.2fx", (float)thread1Stopwatch.getTimeNs()/tasks2Stopwatch.getTimeNs()));
		
		tasks.stopAndWait(10000);
		tasks.start(4);
		
		System.out.println("\nBenchmarking Emat calculation, 4 task threads...");
		Stopwatch tasks4Stopwatch = new Stopwatch().start();
		EnergyMatrix emat4 = new SimpleEnergyMatrixCalculator(ecalc).calcEnergyMatrix(tasks);
		checkEmat(search.emat, emat4);
		tasks4Stopwatch.stop();
		System.out.println(String.format("Speedup: %.2fx", (float)thread1Stopwatch.getTimeNs()/tasks4Stopwatch.getTimeNs()));
		
		tasks.stopAndWait(10000);
	}

	private static void checkEmat(EnergyMatrix exp, EnergyMatrix obs) {
		
		// minimized energies aren't super precise
		final double Epsilon = 1e-3;
		
		for (int pos1=0; pos1<exp.getNumPos(); pos1++) {
			for (int rc1=0; rc1<exp.getNumConfAtPos(pos1); rc1++) {
			
				// singles
				assertThat(obs.getOneBody(pos1, rc1), isRelatively(exp.getOneBody(pos1, rc1), Epsilon));
				
				// pairs
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<exp.getNumConfAtPos(pos2); rc2++) {
						
						assertThat(obs.getPairwise(pos1, rc1, pos2, rc2), isRelatively(exp.getPairwise(pos1, rc1, pos2, rc2), Epsilon));
					}
				}
			}
		}
	}
}
