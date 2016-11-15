package edu.duke.cs.osprey.minimization;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.ForcefieldInteractionsGenerator;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tupexp.LUTESettings;

@SuppressWarnings("unused")
public class BenchmarkMinimization extends TestBase {
	
	public static void main(String[] args)
	throws Exception {
		
		initDefaultEnvironment();
		
		// for these small problems, more than one thread is actually slower
		MultiTermEnergyFunction.setNumThreads(1);
		
		// try to minimize context switches that can upset our timings
		Thread.currentThread().setPriority(Thread.MAX_PRIORITY);
		
		// make a search problem
		System.out.println("Building search problem...");
		
		ResidueFlexibility resFlex = new ResidueFlexibility();
		//resFlex.addMutable("39 43 47 52", "ALA GLY VAL ILE TRP HIS");
		resFlex.addMutable("39 43", "ALA");
		resFlex.addFlexible("40 41 42 44 45");
		boolean doMinimize = true;
		boolean addWt = false;
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
			resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null,
			false, new ArrayList<>()
		);
		
		// calc the energy matrix
		File ematFile = new File("/tmp/benchmarkMinimization.emat.dat");
		if (ematFile.exists()) {
			search.emat = (EnergyMatrix)ObjectIO.readObject(ematFile.getAbsolutePath(), false);
		} else {
			ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
			tasks.start(2);
			SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(EnvironmentVars.curEFcnGenerator, search.confSpace, search.shellResidues);
			search.emat = new SimpleEnergyMatrixCalculator(ecalc).calcEnergyMatrix(tasks);
			tasks.stop();
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
		}
		
		// settings
		final int numConfs = 16;//64;//256;
		
		// get a few arbitrary conformations
		System.out.println("getting confs...");
		search.pruneMat = new PruningMatrix(search.confSpace, 1000);
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 4, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		List<ScoredConf> confs = new ArrayList<>();
		for (int i=0; i<numConfs; i++) {
			confs.add(tree.nextConf());
		}
		
		/* TEMP: sometimes, we just want the ith conf
		// but be sure to ignore the energy warnings, because the order won't match anymore
		List<ScoredConf> newConfs = new ArrayList<>();
		for (int i=0; i<8; i++) {
			newConfs.add(confs.get(0));
		}
		confs = newConfs;
		*/
		
		System.out.println("benchmarking...");
		
		//benchmarkSerial(search, confs);
		benchmarkParallel(search, confs);
	}
	
	private static class Factories {
		
		public final GpuQueuePool openclPool;
		public final GpuStreamPool cudaPool;
		public final Factory<Minimizer,MoleculeModifierAndScorer> originalMinimizers;
		public final Factory<Minimizer,MoleculeModifierAndScorer> simpleMinimizers;
		public final Factory<Minimizer,MoleculeModifierAndScorer> cudaMinimizers;
		public final Factory<EnergyFunction,Molecule> cpuEfuncs;
		public final Factory<EnergyFunction,Molecule> openclEfuncs;
		public final Factory<EnergyFunction,Molecule> cudaEfuncs;
		public final Factory<EnergyFunction,Molecule> bigEfuncs;
		
		public Factories(SearchProblem search, int maxNumStreams) {
			
			System.out.println("making gpu queues...");
			
			openclPool = new GpuQueuePool(1, maxNumStreams);
			
			System.out.println("making gpu streams...");
			
			cudaPool = new GpuStreamPool(1, maxNumStreams);
			
			System.out.println("making minimizers...");
			
			// make minimizer factories
			originalMinimizers = new Factory<Minimizer,MoleculeModifierAndScorer>() {
				@Override
				public Minimizer make(MoleculeModifierAndScorer mof) {
					return new CCDMinimizer(mof, true);
				}
			};
			simpleMinimizers = new Factory<Minimizer,MoleculeModifierAndScorer>() {
				@Override
				public Minimizer make(MoleculeModifierAndScorer mof) {
					return new SimpleCCDMinimizer(mof, new Factory<LineSearcher,Void>() {
						@Override
						public LineSearcher make(Void ignore) {
							return new SurfingLineSearcher();
						}
					});
				}
			};
			
			cudaMinimizers = new Factory<Minimizer,MoleculeModifierAndScorer>() {
				@Override
				public Minimizer make(MoleculeModifierAndScorer mof) {
					return new CudaCCDMinimizer(cudaPool, mof);
				}
			};
			
			System.out.println("making efuncs...");
			
			// make efuncs
			cpuEfuncs = new Factory<EnergyFunction,Molecule>() {
				@Override
				public EnergyFunction make(Molecule mol) {
					return EnvironmentVars.curEFcnGenerator.fullConfEnergy(search.confSpace, search.shellResidues, mol);
				}
			};
			
			GpuEnergyFunctionGenerator openclEgen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), openclPool);
			openclEfuncs = new Factory<EnergyFunction,Molecule>() {
				@Override
				public EnergyFunction make(Molecule mol) {
					return openclEgen.fullConfEnergy(search.confSpace, search.shellResidues, mol);
				}
			};
			
			GpuEnergyFunctionGenerator cudaEgen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), cudaPool);
			cudaEfuncs = new Factory<EnergyFunction,Molecule>() {
				@Override
				public EnergyFunction make(Molecule mol) {
					return cudaEgen.fullConfEnergy(search.confSpace, search.shellResidues, mol);
				}
			};
			
			ForcefieldInteractionsGenerator ffintergen = new ForcefieldInteractionsGenerator();
			bigEfuncs = new Factory<EnergyFunction,Molecule>() {
				@Override
				public EnergyFunction make(Molecule mol) {
					ForcefieldInteractions interactions = ffintergen.makeFullConf(search.confSpace, search.shellResidues, mol);
					return new BigForcefieldEnergy(EnvironmentVars.curEFcnGenerator.ffParams, interactions, BufferTools.Type.Direct);
				}
			};
			
			System.out.println("factories complete");
		}
		
		public void cleanup() {
			openclPool.cleanup();
			cudaPool.cleanup();
		}
	}
	
	private static void benchmarkSerial(SearchProblem search, List<ScoredConf> confs)
	throws Exception {

		Factories f = new Factories(search, 1);
		List<EnergiedConf> minimizedConfs;
		
		System.out.println("\nbenchmarking CPU original...");
		Stopwatch cpuOriginalStopwatch = new Stopwatch().start();
		minimizedConfs = new ConfMinimizer(f.originalMinimizers).minimize(confs, f.cpuEfuncs, search.confSpace);
		System.out.println("precise timing: " + cpuOriginalStopwatch.stop().getTime(TimeUnit.MILLISECONDS));
		checkEnergies(minimizedConfs);
		
		System.out.println("\nbenchmarking CPU simple...");
		Stopwatch cpuSimpleStopwatch = new Stopwatch().start();
		minimizedConfs = new ConfMinimizer(f.simpleMinimizers).minimize(confs, f.cpuEfuncs, search.confSpace);
		System.out.print("precise timing: " + cpuSimpleStopwatch.stop().getTime(TimeUnit.MILLISECONDS));
		System.out.println(String.format(", speedup: %.2fx",
			(double)cpuOriginalStopwatch.getTimeNs()/cpuSimpleStopwatch.getTimeNs()
		));
		checkEnergies(minimizedConfs);
		
		System.out.println("\nbenchmarking CPU simple big...");
		Stopwatch cpuSimpleBigStopwatch = new Stopwatch().start();
		minimizedConfs = new ConfMinimizer(f.simpleMinimizers).minimize(confs, f.bigEfuncs, search.confSpace);
		System.out.print("precise timing: " + cpuSimpleBigStopwatch.stop().getTime(TimeUnit.MILLISECONDS));
		System.out.println(String.format(", speedup: %.2fx",
			(double)cpuOriginalStopwatch.getTimeNs()/cpuSimpleBigStopwatch.getTimeNs()
		));
		checkEnergies(minimizedConfs);
		
		{
			// warm up the energy function
			// avoids anomalous timings on the first conf
			GpuForcefieldEnergy efunc = (GpuForcefieldEnergy)f.openclEfuncs.make(search.confSpace.m);
			efunc.getEnergy();
			efunc.cleanup();
		}
		
		System.out.println("\nbenchmarking OpenCL simple...");
		Stopwatch openclSimpleStopwatch = new Stopwatch().start();
		minimizedConfs = new ConfMinimizer(f.simpleMinimizers).minimize(confs, f.openclEfuncs, search.confSpace);
		System.out.print("precise timing: " + openclSimpleStopwatch.stop().getTime(TimeUnit.MILLISECONDS));
		System.out.println(String.format(", speedup over CPU simple: %.2fx",
			(double)cpuSimpleStopwatch.getTimeNs()/openclSimpleStopwatch.getTimeNs()
		));
		checkEnergies(minimizedConfs);
		
		{
			// warm up the energy function
			// avoids anomalous timings on the first conf
			GpuForcefieldEnergy efunc = (GpuForcefieldEnergy)f.cudaEfuncs.make(search.confSpace.m);
			efunc.getEnergy();
			efunc.cleanup();
		}
		
		System.out.println("\nbenchmarking Cuda simple...");
		Stopwatch cudaSimpleStopwatch = new Stopwatch().start();
		minimizedConfs = new ConfMinimizer(f.simpleMinimizers).minimize(confs, f.cudaEfuncs, search.confSpace);
		System.out.print("precise timing: " + cudaSimpleStopwatch.stop().getTime(TimeUnit.MILLISECONDS));
		System.out.println(String.format(", speedup over CPU simple: %.2fx, speedup over OpenCL simple: %.2fx",
			(double)cpuSimpleStopwatch.getTimeNs()/cudaSimpleStopwatch.getTimeNs(),
			(double)openclSimpleStopwatch.getTimeNs()/cudaSimpleStopwatch.getTimeNs()
		));
		checkEnergies(minimizedConfs);
		
		// TODO: warm up CCD?
		
		System.out.println("\nbenchmarking Cuda CCD..");
		Stopwatch cudaCCDStopwatch = new Stopwatch().start();
		minimizedConfs = new ConfMinimizer(f.cudaMinimizers).minimize(confs, f.bigEfuncs, search.confSpace);
		System.out.print("precise timing: " + cudaCCDStopwatch.stop().getTime(TimeUnit.MILLISECONDS));
		System.out.println(String.format(", speedup over CPU simple: %.2fx, speedup over Cuda simple: %.2fx",
			(double)cpuSimpleStopwatch.getTimeNs()/cudaCCDStopwatch.getTimeNs(),
			(double)cudaSimpleStopwatch.getTimeNs()/cudaCCDStopwatch.getTimeNs()
		));
		checkEnergies(minimizedConfs);
		
		// cleanup
		f.cleanup();
	}
	
	private static void benchmarkParallel(SearchProblem search, List<ScoredConf> confs)
	throws Exception {
		
		// settings
		final int[] numThreadsList = { 1 };//, 2, 4 };
		final int[] numStreamsList = { 1, 2, 4, 8 };//, 16, 32 };
		int maxNumStreams = numStreamsList[numStreamsList.length - 1];
		
		Factories f = new Factories(search, maxNumStreams);
		List<EnergiedConf> minimizedConfs;
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		
		// benchmark cpu
		Stopwatch oneCpuStopwatch = null;
		for (int numThreads : numThreadsList) {
			
			tasks.start(numThreads);
			
			System.out.println("\nBenchmarking " + numThreads + " thread(s) with CPU efuncs...");
			Stopwatch taskStopwatch = new Stopwatch().start();
			minimizedConfs = new ConfMinimizer(f.simpleMinimizers).minimize(confs, f.cpuEfuncs, search.confSpace, tasks);
			taskStopwatch.stop();
			if (oneCpuStopwatch == null) {
				oneCpuStopwatch = taskStopwatch;
			}
			System.out.println("precise timing: " + taskStopwatch.getTime(TimeUnit.MILLISECONDS));
			System.out.println(String.format("Speedup over one: %.2fx",
				(float)oneCpuStopwatch.getTimeNs()/taskStopwatch.getTimeNs()
			));
			checkEnergies(minimizedConfs);
			
			tasks.stopAndWait(10000);
		}
		
		// benchmark opencl
		Stopwatch oneOpenCLStopwatch = null;
		for (int numThreads : numThreadsList) {
			
			tasks.start(numThreads);
			
			System.out.println("\nBenchmarking " + numThreads + " thread(s) with OpenCL efuncs...");
			Stopwatch taskStopwatch = new Stopwatch().start();
			minimizedConfs = new ConfMinimizer(f.simpleMinimizers).minimize(confs, f.openclEfuncs, search.confSpace, tasks);
			taskStopwatch.stop();
			if (oneOpenCLStopwatch == null) {
				oneOpenCLStopwatch = taskStopwatch;
			}
			System.out.println("precise timing: " + taskStopwatch.getTime(TimeUnit.MILLISECONDS));
			System.out.println(String.format("Speedup over one CPU: %.2fx, Speedup over one: %.2fx",
				(float)oneCpuStopwatch.getTimeNs()/taskStopwatch.getTimeNs(),
				(float)oneOpenCLStopwatch.getTimeNs()/taskStopwatch.getTimeNs()
			));
			checkEnergies(minimizedConfs);
			
			tasks.stopAndWait(10000);
		}
		
		// benchmark cuda
		Stopwatch oneCudaStopwatch = null;
		for (int numThreads : numThreadsList) {
			
			tasks.start(numThreads);
			
			System.out.println("\nBenchmarking " + numThreads + " thread(s) with Cuda efuncs...");
			Stopwatch taskStopwatch = new Stopwatch().start();
			minimizedConfs = new ConfMinimizer(f.simpleMinimizers).minimize(confs, f.cudaEfuncs, search.confSpace, tasks);
			taskStopwatch.stop();
			if (oneCudaStopwatch == null) {
				oneCudaStopwatch = taskStopwatch;
			}
			System.out.println("precise timing: " + taskStopwatch.getTime(TimeUnit.MILLISECONDS));
			System.out.println(String.format("Speedup over one CPU: %.2fx, Speedup over one: %.2fx",
				(float)oneCpuStopwatch.getTimeNs()/taskStopwatch.getTimeNs(),
				(float)oneCudaStopwatch.getTimeNs()/taskStopwatch.getTimeNs()
			));
			checkEnergies(minimizedConfs);
			
			tasks.stopAndWait(10000);
		}

		// benchmark cuda ccd
		Stopwatch oneCudaCCDStopwatch = null;
		for (int numStreams : numStreamsList) {
			
			tasks.start(numStreams);
			
			System.out.println("\nBenchmarking " + numStreams + " CCD stream(s)...");
			Stopwatch taskStopwatch = new Stopwatch().start();
			minimizedConfs = new ConfMinimizer(f.cudaMinimizers).minimize(confs, f.bigEfuncs, search.confSpace, tasks);
			taskStopwatch.stop();
			if (oneCudaCCDStopwatch == null) {
				oneCudaCCDStopwatch = taskStopwatch;
			}
			System.out.println("precise timing: " + taskStopwatch.getTime(TimeUnit.MILLISECONDS));
			System.out.println(String.format("Speedup over one CPU: %.2fx, Speedup over one: %.2fx",
				(float)oneCpuStopwatch.getTimeNs()/taskStopwatch.getTimeNs(),
				(float)oneCudaCCDStopwatch.getTimeNs()/taskStopwatch.getTimeNs()
			));
			checkEnergies(minimizedConfs);
			
			tasks.stopAndWait(10000);
		}
		
		f.cleanup();
	}

	private static void checkEnergies(List<EnergiedConf> minimizedConfs) {
		
		// what do we expect the energies to be?
		final double[] expectedEnergies = {
		  -89.40966965,   -89.10792020,   -89.80959784,   -88.63999170,   -89.12813427,   -89.50404295,   -88.39619839,   -88.88944800,
		  -88.91538573,   -88.37400350,   -88.72521732,   -88.95813013,   -88.71957443,   -89.13542393,   -88.39342982,   -88.61512253,
		  -87.67170876,   -87.75762493,   -88.82437753,   -87.49111955,   -88.39267148,   -88.73093938,   -88.16624474,   -88.47911298,
		  -88.77040379,   -88.59812445,   -88.63729374,   -88.83029528,   -88.07022384,   -87.82115594,   -89.15034617,   -89.52776262,
		  -88.08922551,   -87.24538989,   -87.66296341,   -89.47261229,   -88.53548693,   -88.21416862,   -87.18239056,   -88.37126489,
		  -88.91533055,   -88.95432276,   -88.34024189,   -88.53617041,   -87.76065876,   -87.75246826,   -89.32887293,   -89.12214183,
		  -87.53435849,   -89.30674536,   -88.86121108,   -88.00498514,   -89.24745408,   -86.93536186,   -87.83485265,   -89.18378421,
		  -87.60530136,   -87.88059458,   -88.99239407,   -89.00570101,   -88.47514883,   -88.62549053,   -89.05482774,   -88.65430730
		};
		
		final double Epsilon = 1e-3;
		
		int n = minimizedConfs.size();
		for (int i=0; i<n; i++) {
			
			double energy = minimizedConfs.get(i).getEnergy();
			
			if (i < expectedEnergies.length) {
				
				if (Double.isNaN(energy)) {
					System.out.println(String.format("\tWARNING: invalid energy: i:%-3d  exp:%12.8f  obs: NaN",
						i, expectedEnergies[i]
					));
				} else {
				
					double absErr = energy - expectedEnergies[i];
					if (absErr > Epsilon) {
						
						System.out.println(String.format("\tWARNING: low precision energy: i:%-3d  exp:%12.8f  obs:%12.8f       absErr:%12.8f",
							i, expectedEnergies[i], energy, absErr
						));
						
					} else if (absErr < -Epsilon) {
					
						System.out.println(String.format("\t              improved energy: i:%-3d  exp:%12.8f  obs:%12.8f  improvement:%12.8f",
							i, expectedEnergies[i], energy, -absErr
						));
					}
				}
				
			} else {
				
				// print the energy so we can write the accuracy test
				if (i > 0) {
					System.out.print(",");
				}
				System.out.print(i % 8 == 0 ? "\n" : " ");
				System.out.print(String.format("%14.8f", energy));
			}
		}
	}
}
