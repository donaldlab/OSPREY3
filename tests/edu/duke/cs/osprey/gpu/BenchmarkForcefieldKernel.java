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
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.TimeFormatter;
import edu.duke.cs.osprey.tools.TimingThread;

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
		//String mutRes = "39 43";
		String mutRes = "39 43 46 47";
		//String flexRes = "";
		//String flexRes = "40";
		//String flexRes = "40 41";
		//String flexRes = "40 41 42 44 45";
		String flexRes = "40 41 42 44 45 48 49 50 51 52 53";
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
		GpuEnergyFunctionGenerator gpuegen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), new GpuQueuePool(8, 2));
		
		benchmarkEfunc(search, egen, gpuegen);
		//benchmarkEmat(search, egen, gpuegen);
		//benchmarkMinimize(search, egen, gpuegen);
	}
	
	private static void benchmarkEfunc(SearchProblem search, EnergyFunctionGenerator egen, GpuEnergyFunctionGenerator gpuegen)
	throws Exception {
		
		// figure out how many threads to use based on the gpu queue pool
		List<Integer> numThreadsList = new ArrayList<>();
		for (int i=1; i<=gpuegen.getQueuePool().getNumQueues(); i<<=1) {
			numThreadsList.add(i);
		}
		
		System.out.println("\nFull conf energy:");
		benchmarkEfunc(
			//500,
			//1000,
			//2000,
			10000,
			search.confSpace.m,
			new Factory<EnergyFunction,Molecule>() {
				@Override
				public EnergyFunction make(Molecule mol) {
					return egen.fullConfEnergy(search.confSpace, search.shellResidues, mol);
				}
			},
			new Factory<GpuForcefieldEnergy,Molecule>() {
				@Override
				public GpuForcefieldEnergy make(Molecule mol) {
					return gpuegen.fullConfEnergy(search.confSpace, search.shellResidues, mol);
				}
			},
			numThreadsList
		);
		
		/* TEMP
		System.out.println("\nIntra and shell energy:");
		benchmarkEfunc(
			10000,
			search.confSpace.m,
			new Factory<EnergyFunction,Molecule>() {
				@Override
				public EnergyFunction make(Molecule mol) {
					return egen.intraAndShellEnergy(search.confSpace.posFlex.get(0).res, search.shellResidues, mol); 
				}
			},
			new Factory<GpuForcefieldEnergy,Molecule>() {
				@Override
				public GpuForcefieldEnergy make(Molecule mol) {
					return gpuegen.intraAndShellEnergy(search.confSpace.posFlex.get(0).res, search.shellResidues, mol);
				}
			},
			numThreadsList
		);
		
		System.out.println("\nPairwise energy:");
		// TODO: GPU is actually significantly slower for these terms
		// there's so few atom pairs, the overhead with the GPU is slowing us down
		// need to optimize more, maybe look into faster memory transfers?
		// NOTE: profiling says the memory transfers are really fast, ~4/43 us or ~9%
		// most of the overhead seems to be coming from synchronization with the GPU, ~26/43 us or ~60%
		// don't think there's anything we can do to speed that up...
		// sync overhead is relatively smaller for other sizes, ~18% for full conf energy, ~42% for intra and shell energy
		benchmarkEfunc(
			50000,
			search.confSpace.m,
			new Factory<EnergyFunction,Molecule>() {
				@Override
				public EnergyFunction make(Molecule mol) {
					return egen.resPairEnergy(search.confSpace.posFlex.get(0).res, search.confSpace.posFlex.get(2).res, mol); 
				}
			},
			new Factory<GpuForcefieldEnergy,Molecule>() {
				@Override
				public GpuForcefieldEnergy make(Molecule mol) {
					return gpuegen.resPairEnergy(search.confSpace.posFlex.get(0).res, search.confSpace.posFlex.get(2).res, mol);
				}
			},
			numThreadsList
		);
		*/
	}
	
	private static void benchmarkEfunc(int numRuns, Molecule baseMol, Factory<EnergyFunction,Molecule> efuncs, Factory<GpuForcefieldEnergy,Molecule> gpuefuncs, List<Integer> numThreadsList) {
		
		// benchmark cpus
		Stopwatch cpuStopwatch = null;
		for (int numThreads : numThreadsList) {
			System.out.print("Benchmarking " + numThreads + " CPUs... ");
			
			// make the thread pool
			List<TimingThread> threads = new ArrayList<>();
			for (int i=0; i<numThreads; i++) {
				threads.add(new TimingThread("Cpu-" + i) {
					
					private Molecule mol;
					private EnergyFunction efunc;
					
					@Override
					public void warmup() {
						mol = new Molecule(baseMol);
						efunc = efuncs.make(mol);
						for (int k=0; k<10; k++) {
							efunc.getEnergy();
						}
					}
					
					@Override
					public void time() {
						int numLocalRuns = numRuns/numThreads;
						for (int k=0; k<numLocalRuns; k++) {
							efunc.getEnergy();
						}
					}
				});
			}
			
			// run threads and wait
			Stopwatch cpusStopwatch = TimingThread.timeThreads(threads);
			
			if (numThreads == 1) {
				cpuStopwatch = cpusStopwatch;
			}
			
			// show timing info
			System.out.println(String.format(" finished in %s, avg time per op: %s, speedup: %.2fx",
				cpusStopwatch.getTime(2),
				TimeFormatter.format(cpusStopwatch.getTimeNs()/numRuns, TimeUnit.MICROSECONDS),
				(double)cpuStopwatch.getTimeNs()/cpusStopwatch.getTimeNs()
			));
		}
		
		// benchmark gpus
		for (int numThreads : numThreadsList) {
			System.out.print("Benchmarking " + numThreads + " GPUs... ");
			
			// TEMP
			//System.out.println();
			
			// make the thread pool
			List<TimingThread> threads = new ArrayList<>();
			for (int i=0; i<numThreads; i++) {
				threads.add(new TimingThread("Gpu-" + i) {
					
					private Molecule mol;
					private GpuForcefieldEnergy gpuefunc;
					
					@Override
					public void warmup() {
						mol = new Molecule(baseMol);
						gpuefunc = gpuefuncs.make(mol);
						for (int j=0; j<10; j++) {
							gpuefunc.getEnergy();
						}
						/* TEMP
						System.out.println(String.format("Thread %s has gpu %d",
							getName(), efunc.getKernel().getQueue().getDevice().getID()
						));
						*/
					}
					
					@Override
					public void time() {
						
						int numLocalRuns = numRuns/numThreads;
						for (int j=0; j<numLocalRuns; j++) {
							gpuefunc.getEnergy();
						}
						
						gpuefunc.cleanup();
					}
				});
			}
			
			// run the threads and wait
			Stopwatch gpuStopwatch = TimingThread.timeThreads(threads);
			
			// show timing info
			System.out.println(String.format(" finished in %s, avg time per op: %s, speedup: %.2fx",
				gpuStopwatch.getTime(2),
				TimeFormatter.format(gpuStopwatch.getTimeNs()/numRuns, TimeUnit.MICROSECONDS),
				(double)cpuStopwatch.getTimeNs()/gpuStopwatch.getTimeNs()
			));
		}
		
		// do a final gpu run for profiling
		GpuForcefieldEnergy gpuefunc = gpuefuncs.make(null);
		if (gpuefunc.getKernel().getQueue().isProfilingEnabled()) {
			gpuefunc.startProfile();
		}
		gpuefunc.getEnergy();
		System.out.println("GPU profiling info:");
		System.out.println("atom pairs:      " + gpuefunc.getForcefieldEnergy().getNumAtomPairs());
		System.out.println("GPU memory used: " + gpuefunc.getKernel().getGpuBytesNeeded()/1024 + " KiB");
		if (gpuefunc.getKernel().getQueue().isProfilingEnabled()) {
			System.out.print(gpuefunc.dumpProfile());
		}
		gpuefunc.cleanup();
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
