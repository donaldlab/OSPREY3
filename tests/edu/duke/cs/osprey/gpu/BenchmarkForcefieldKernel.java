package edu.duke.cs.osprey.gpu;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.TimeFormatter;

public class BenchmarkForcefieldKernel extends TestBase {
	
	// NOTE: useful info for optimizing gpu kernels:
	// http://www.nvidia.com/content/GTC/documents/1068_GTC09.pdf
	
	public static void main(String[] args)
	throws Exception {
		
		initDefaultEnvironment();
		
		// for these small problems, more than one thread is actually slower
		MultiTermEnergyFunction.setNumThreads(1);
		
		// make a search problem
		System.out.println("Building search problem...");
		
		//String aaNames = "ALA VAL LEU ILE PHE TYR TRP CYS MET SER THR LYS ARG HIE HID ASP GLU ASN GLN GLY";
		//String aaNames = "ALA VAL LEU ILE";
		String aaNames = "ALA";
		//String mutRes = "39";
		String mutRes = "39 43";
		//String flexRes = "";
		//String flexRes = "40";
		//String flexRes = "40 41";
		String flexRes = "40 41 42 44 45";
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
		GpuEnergyFunctionGenerator gpuegen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), new GpuQueuePool(2, 2));
		
		//benchmarkEfunc(search, egen, gpuegen);
		//benchmarkEmat(search, egen, gpuegen);
		benchmarkMinimize(search, egen, gpuegen);
	}
	
	private static void benchmarkEfunc(SearchProblem search, EnergyFunctionGenerator egen, GpuEnergyFunctionGenerator gpuegen)
	throws Exception {
		
		List<GpuForcefieldEnergy> gpuefuncs = null;
		
		System.out.println("\nFull conf energy:");
		gpuefuncs = new ArrayList<>();
		for (int i=0; i<gpuegen.getQueuePool().getNumQueues(); i++) {
			gpuefuncs.add(gpuegen.fullConfEnergy(search.confSpace, search.shellResidues));
		}
		benchmarkEfunc(4000,
			egen.fullConfEnergy(search.confSpace, search.shellResidues),
			gpuefuncs
		);
		
		System.out.println("\nIntra and shell energy:");
		gpuefuncs = new ArrayList<>();
		for (int i=0; i<gpuegen.getQueuePool().getNumQueues(); i++) {
			gpuefuncs.add(gpuegen.intraAndShellEnergy(search.confSpace.posFlex.get(0).res, search.shellResidues));
		}
		benchmarkEfunc(20000,
			egen.intraAndShellEnergy(search.confSpace.posFlex.get(0).res, search.shellResidues),
			gpuefuncs
		);
		
		System.out.println("\nPairwise energy:");
		// TODO: GPU is actually significantly slower for these terms
		// there's so few atom pairs, the overhead with the GPU is slowing us down
		// need to optimize more, maybe look into faster memory transfers?
		// NOTE: profiling says the memory transfers are really fast, ~4/43 us or ~9%
		// most of the overhead seems to be coming from synchronization with the GPU, ~26/43 us or ~60%
		// don't think there's anything we can do to speed that up...
		// sync overhead is relatively smaller for other sizes, ~18% for full conf energy, ~42% for intra and shell energy
		gpuefuncs = new ArrayList<>();
		for (int i=0; i<gpuegen.getQueuePool().getNumQueues(); i++) {
			gpuefuncs.add(gpuegen.resPairEnergy(search.confSpace.posFlex.get(0).res, search.confSpace.posFlex.get(2).res));
		}
		benchmarkEfunc(100000,
			egen.resPairEnergy(search.confSpace.posFlex.get(0).res, search.confSpace.posFlex.get(2).res),
			gpuefuncs
		);
	}
	
	private static void benchmarkEfunc(int numRuns, EnergyFunction efunc, List<GpuForcefieldEnergy> gpuefuncs) {
		
		// benchmark the cpu
		System.out.print("Benchmarking CPU...");
		Stopwatch cpuStopwatch = new Stopwatch().start(); 
		for (int i=0; i<numRuns; i++) {
			efunc.getEnergy();
		}
		cpuStopwatch.stop();
		System.out.println(String.format(" finished in %s, avg time per op: %s",
			cpuStopwatch.getTime(2),
			TimeFormatter.format(cpuStopwatch.getTimeNs()/numRuns, TimeUnit.MICROSECONDS)
		));
		
		// benchmark the gpu
		System.out.print("Benchmarking GPU...\n");
		
		// set up thread pool to match queue pool
		List<Thread> threads = new ArrayList<>();
		final List<String> profiles = new ArrayList<>();
		for (int i=0; i<gpuefuncs.size(); i++) {
			final GpuForcefieldEnergy gpuefunc = gpuefuncs.get(i);
			final boolean isFirstThread = i == 0;
			Thread thread = new Thread("Gpu-" + i) {
				@Override
				public void run() {
					
					int numLocalRuns = numRuns/gpuefuncs.size();
					boolean useProfiling = gpuefunc.getKernel().getQueue().isProfilingEnabled();
					
					for (int j=0; j<numLocalRuns; j++) {
						
						boolean isProfileRun = useProfiling && isFirstThread && (j == 0 || j == numLocalRuns - 1);
						if (isProfileRun) {
							gpuefunc.startProfile();
						}
						
						gpuefunc.getEnergy();
						
						if (isProfileRun) {
							profiles.add(gpuefunc.dumpProfile());
						}
					}
				}
			};
			threads.add(thread);
		}
		
		Stopwatch gpuStopwatch = new Stopwatch().start();
		for (Thread thread : threads) {
			thread.start();
		}
		for (Thread thread : threads) {
			try {
				thread.join();
			} catch (InterruptedException ex) {
				throw new Error(ex);
			}
		}
		
		gpuStopwatch.stop();
		System.out.println(String.format(" finished in %s, avg time per op: %s, speedup: %.2fx, numPairs: %d, GPU mem used: %.2f MiB",
			gpuStopwatch.getTime(2),
			TimeFormatter.format(gpuStopwatch.getTimeNs()/numRuns, TimeUnit.MICROSECONDS),
			(double)cpuStopwatch.getTimeNs()/gpuStopwatch.getTimeNs(),
			gpuefuncs.get(0).getForcefieldEnergy().getNumAtomPairs(),
			(double)gpuefuncs.get(0).getKernel().getGpuBytesNeeded()/1024/1024
		));
		if (!profiles.isEmpty()) {
			System.out.println("GPU profiling info:");
			for (int i=0; i<profiles.size(); i++) {
				System.out.print(String.format("%s run:\n\t%s\n", i == 0 ? "first" : "last", profiles.get(i).replace("\n", "\n\t").trim()));
			}
		}
		
		// cleanup
		for (GpuForcefieldEnergy gpuefunc : gpuefuncs) {
			gpuefunc.cleanup();
		}
	}
	
	private static void benchmarkEmat(SearchProblem search, EnergyFunctionGenerator egen, GpuEnergyFunctionGenerator gpuegen)
	throws Exception {
		
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, search.confSpace, search.shellResidues);
		SimpleEnergyCalculator gpuecalc = new SimpleEnergyCalculator(gpuegen, search.confSpace, search.shellResidues);
		
		// benchmark the cpu
		System.out.println("\nBenchmarking CPU...");
		Stopwatch cpuStopwatch = new Stopwatch().start();
		EnergyMatrix emat = new SimpleEnergyMatrixCalculator(ecalc).calcEnergyMatrix();
		cpuStopwatch.stop();
		
		// benchmark the gpu
		System.out.println("\nBenchmarking GPU...");
		Stopwatch gpuStopwatch = new Stopwatch().start();
		EnergyMatrix gpuemat = new SimpleEnergyMatrixCalculator(gpuecalc).calcEnergyMatrix();
		gpuStopwatch.stop();
		
		// calculate speedup
		System.out.println(String.format("\nspeedup: %.2fx", (double)cpuStopwatch.getTimeNs()/gpuStopwatch.getTimeNs()));
		
		// check the result
		// TODO: apparently these results don't match, need to find out why
		for (int pos1=0; pos1<emat.getNumPos(); pos1++) {
			for (int rc1=0; rc1<emat.getNumConfAtPos(pos1); rc1++) {
				
				assertThat(emat.getOneBody(pos1, rc1), isRelatively(gpuemat.getOneBody(pos1, rc1)));
				
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<emat.getNumConfAtPos(pos2); rc2++) {
						
						assertThat(emat.getPairwise(pos1, rc1, pos2, rc2), isRelatively(gpuemat.getPairwise(pos1, rc1, pos2, rc2)));
					}
				}
			}
		}
	}
	
	private static void benchmarkMinimize(SearchProblem search, EnergyFunctionGenerator egen, GpuEnergyFunctionGenerator gpuegen)
	throws Exception {
		
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, search.confSpace, search.shellResidues);
		
		int numConfs = 10;
		
		// get a few arbitrary conformations
		search.emat = new SimpleEnergyMatrixCalculator(ecalc).calcEnergyMatrix();
		search.pruneMat = new PruningMatrix(search.confSpace, 1000);
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 4, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		List<ScoredConf> confs = new ArrayList<>();
		for (int i=0; i<numConfs; i++) {
			confs.add(tree.nextConf());
		}
		
		double energy;
		
		// benchmark the cpu
		System.out.println("\nBenchmarking CPU...");
		EnergyFunction efunc = egen.fullConfEnergy(search.confSpace, search.shellResidues);
		Stopwatch cpuStopwatch = new Stopwatch().start();
		energy = 0;
		for (ScoredConf conf : confs) {
			energy += minimize(efunc, search.confSpace, conf);
		}
		cpuStopwatch.stop();
		System.out.println("\tfinished in " + cpuStopwatch.getTime(2));
		System.out.println("\tenergy sum: " + energy);
		
		// benchmark the gpu
		System.out.println("\nBenchmarking GPU...");
		GpuForcefieldEnergy gpuefunc = gpuegen.fullConfEnergy(search.confSpace, search.shellResidues);
		Stopwatch gpuStopwatch = new Stopwatch().start();
		energy = 0;
		for (ScoredConf conf : confs) {
			energy += minimize(gpuefunc, search.confSpace, conf);
		}
		gpuStopwatch.stop();
		gpuefunc.cleanup();
		System.out.println("\tfinished in " + gpuStopwatch.getTime(2));
		System.out.println("\tenergy: " + energy);
		
		// calculate speedup
		System.out.println(String.format("\nspeedup: %.2fx", (double)cpuStopwatch.getTimeNs()/gpuStopwatch.getTimeNs())); 
	}
	
	private static double minimize(EnergyFunction efunc, ConfSpace confSpace, ScoredConf conf) {
		RCTuple tuple = new RCTuple(conf.getAssignments());
		MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, confSpace, tuple);
		new CCDMinimizer(mof, true).minimize();
		//new SimpleCCDMinimizer(mof).minimize();
		return efunc.getEnergy();
	}
}
