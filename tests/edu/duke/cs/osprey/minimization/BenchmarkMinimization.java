package edu.duke.cs.osprey.minimization;

import java.io.File;
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
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
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
import edu.duke.cs.osprey.gpu.GpuQueuePool;
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
		resFlex.addMutable("39 43", "ALA");
		resFlex.addFlexible("40 41 42 44 45");
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
		final int numConfs = 64;//256;
		
		// get a few arbitrary conformations
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
		
		//benchmarkParallelism(search, confs);
		benchmarkGpu(search, confs);
	}
	
	private static void benchmarkParallelism(SearchProblem search, List<ScoredConf> confs)
	throws Exception {
		
		// settings
		final int[] numThreadsList = { 1, 2, 4, 8 };//, 16 };
		final boolean useGpu = false;
		
		int maxNumThreads = numThreadsList[numThreadsList.length - 1];
		
		// get the energy function generator
		final EnergyFunctionGenerator egen;
		if (useGpu) {
			GpuQueuePool gpuPool = new GpuQueuePool(maxNumThreads, 1);
			//GpuQueuePool gpuPool = new GpuQueuePool(1, maxNumThreads);
			egen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), gpuPool);
		} else {
			egen = EnvironmentVars.curEFcnGenerator;
		}
		
		// make the energy function factory
		Factory<EnergyFunction,Molecule> efuncs = new Factory<EnergyFunction,Molecule>() {
			@Override
			public EnergyFunction make(Molecule mol) {
				return egen.fullConfEnergy(search.confSpace, search.shellResidues, mol);
			}
		};
		
		List<EnergiedConf> minimizedConfs;
		
		// benchmark base minimization
		System.out.println("\nBenchmarking main thread...");
		Stopwatch baseStopwatch = new Stopwatch().start();
		minimizedConfs = new ConfMinimizer().minimize(confs, efuncs, search.confSpace);
		baseStopwatch.stop();
		System.out.println("precise timing: " + baseStopwatch.getTime(TimeUnit.MILLISECONDS));
		checkEnergies(minimizedConfs);
		
		// benchmark parallel minimization
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		
		for (int numThreads : numThreadsList) {
			
			tasks.start(numThreads);
			
			System.out.println("\nBenchmarking " + numThreads + " task thread(s)...");
			Stopwatch taskStopwatch = new Stopwatch().start();
			minimizedConfs = new ConfMinimizer().minimize(confs, efuncs, search.confSpace, tasks);
			taskStopwatch.stop();
			System.out.println("precise timing: " + taskStopwatch.getTime(TimeUnit.MILLISECONDS));
			System.out.println(String.format("Speedup: %.2fx", (float)baseStopwatch.getTimeNs()/taskStopwatch.getTimeNs()));
			checkEnergies(minimizedConfs);
			
			tasks.stopAndWait(10000);
		}
	}

	private static void benchmarkGpu(SearchProblem search, List<ScoredConf> confs)
	throws Exception {
		
		// make minimizer factories
		Factory<Minimizer,MoleculeModifierAndScorer> originalMinimizers = new Factory<Minimizer,MoleculeModifierAndScorer>() {
			@Override
			public Minimizer make(MoleculeModifierAndScorer mof) {
				return new CCDMinimizer(mof, true);
			}
		};
		Factory<Minimizer,MoleculeModifierAndScorer> simpleMinimizers = new Factory<Minimizer,MoleculeModifierAndScorer>() {
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
		Factory<Minimizer,MoleculeModifierAndScorer> pipelinedMinimizers = new Factory<Minimizer,MoleculeModifierAndScorer>() {
			@Override
			public Minimizer make(MoleculeModifierAndScorer mof) {
				return new SimpleCCDMinimizer(mof, new Factory<LineSearcher,Void>() {
					@Override
					public LineSearcher make(Void ignore) {
						return new GpuStyleSurfingLineSearcher();
					}
				});
			}
		};
		Factory<Minimizer,MoleculeModifierAndScorer> gpuPipelinedMinimizers = new Factory<Minimizer,MoleculeModifierAndScorer>() {
			@Override
			public Minimizer make(MoleculeModifierAndScorer mof) {
				return new SimpleCCDMinimizer(mof, new Factory<LineSearcher,Void>() {
					@Override
					public LineSearcher make(Void ignore) {
						return new GpuSurfingLineSearcher();
					}
				});
			}
		};
		
		// make efuncs
		Factory<EnergyFunction,Molecule> efuncs = new Factory<EnergyFunction,Molecule>() {
			@Override
			public EnergyFunction make(Molecule mol) {
				return EnvironmentVars.curEFcnGenerator.fullConfEnergy(search.confSpace, search.shellResidues, mol);
			}
		};
		GpuQueuePool gpuPool = new GpuQueuePool(1);
		GpuEnergyFunctionGenerator gpuegen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), gpuPool);
		Factory<EnergyFunction,Molecule> gpuefuncs = new Factory<EnergyFunction,Molecule>() {
			@Override
			public EnergyFunction make(Molecule mol) {
				return gpuegen.fullConfEnergy(search.confSpace, search.shellResidues, mol);
			}
		};
		
		// warm up the energy function
		// avoids anomolous timings on the first conf
		EnergyFunction efunc = gpuefuncs.make(search.confSpace.m);
		efunc.getEnergy();
		((EnergyFunction.NeedsCleanup)efunc).cleanup();
		
		List<EnergiedConf> minimizedConfs;
		
		System.out.println("\nbenchmarking CPU original...");
		Stopwatch cpuOriginalStopwatch = new Stopwatch().start();
		minimizedConfs = new ConfMinimizer(originalMinimizers).minimize(confs, efuncs, search.confSpace);
		System.out.println("precise timing: " + cpuOriginalStopwatch.stop().getTime(TimeUnit.MILLISECONDS));
		checkEnergies(minimizedConfs);
		
		System.out.println("\nbenchmarking CPU simple...");
		Stopwatch cpuSimpleStopwatch = new Stopwatch().start();
		minimizedConfs = new ConfMinimizer(simpleMinimizers).minimize(confs, efuncs, search.confSpace);
		System.out.print("precise timing: " + cpuSimpleStopwatch.stop().getTime(TimeUnit.MILLISECONDS));
		System.out.println(String.format(", speedup: %.2fx", (double)cpuOriginalStopwatch.getTimeNs()/cpuSimpleStopwatch.getTimeNs()));
		checkEnergies(minimizedConfs);
		
		System.out.println("\nbenchmarking CPU pipelined...");
		Stopwatch cpuPipelinedStopwatch = new Stopwatch().start();
		minimizedConfs = new ConfMinimizer(pipelinedMinimizers).minimize(confs, efuncs, search.confSpace);
		System.out.print("precise timing: " + cpuPipelinedStopwatch.stop().getTime(TimeUnit.MILLISECONDS));
		System.out.println(String.format(", speedup: %.2fx", (double)cpuOriginalStopwatch.getTimeNs()/cpuPipelinedStopwatch.getTimeNs()));
		checkEnergies(minimizedConfs);
		
		System.out.println("\nbenchmarking GPU original...");
		Stopwatch gpuOriginalStopwatch = new Stopwatch().start();
		minimizedConfs = new ConfMinimizer(originalMinimizers).minimize(confs, gpuefuncs, search.confSpace);
		System.out.print("precise timing: " + gpuOriginalStopwatch.stop().getTime(TimeUnit.MILLISECONDS));
		System.out.println(String.format(", speedup: %.2fx", (double)cpuOriginalStopwatch.getTimeNs()/gpuOriginalStopwatch.getTimeNs()));
		checkEnergies(minimizedConfs);
		
		System.out.println("\nbenchmarking GPU pipelined...");
		Stopwatch gpuPipelinedStopwatch = new Stopwatch().start();
		minimizedConfs = new ConfMinimizer(gpuPipelinedMinimizers).minimize(confs, gpuefuncs, search.confSpace);
		System.out.print("precise timing: " + gpuPipelinedStopwatch.stop().getTime(TimeUnit.MILLISECONDS));
		System.out.println(String.format(", speedup over cpu: %.2fx, speedup over original gpu: %.2fx",
			(double)cpuOriginalStopwatch.getTimeNs()/gpuPipelinedStopwatch.getTimeNs(),
			(double)gpuOriginalStopwatch.getTimeNs()/gpuPipelinedStopwatch.getTimeNs()
		));
		checkEnergies(minimizedConfs);
	}
	
	private static void checkEnergies(List<EnergiedConf> minimizedConfs) {
		
		// what do we expect the energies to be?
		final double[] expectedEnergies = {
			 -107.01471459,  -107.14427777,  -106.79145720,  -106.92142505,  -106.28769323,  -107.11806915,  -106.67892491,  -106.41908246,
			 -106.89600270,  -106.74467994,  -106.45688996,  -106.95533586,  -106.43785971,  -106.52194061,  -106.39162552,  -106.72599183,
			 -106.10055035,  -106.73371624,  -106.45589638,  -105.95183007,  -106.40589795,  -106.16492360,  -106.50474688,  -105.73727276,
			 -105.82709238,  -106.01931530,  -106.23489220,  -106.38648271,  -106.23062795,  -106.01634700,  -106.51222044,  -106.13624516,
			 -105.51519164,  -105.58612419,  -106.16496207,  -105.99895146,  -105.72750978,  -105.98690390,  -106.18174766,  -105.79283124,
			 -106.28946688,  -105.72932750,  -106.31993137,  -105.46572215,  -105.37426006,  -105.64006776,  -106.15810861,  -105.00807833,
			 -105.50647942,  -105.76577844,  -105.95677321,  -105.30211004,  -106.10569822,  -106.34665576,  -104.96297178,  -105.24942710,
			 -105.65937038,  -105.80862838,  -106.04288092,  -105.28676434,  -105.78599604,  -105.60879882,  -105.59757051,  -105.93243523,
			 -104.85660318,  -105.82786054,  -105.52743753,  -106.07292822,  -104.76277318,  -105.70197187,  -105.41907447,  -105.59816112,
			 -105.57514850,  -105.00078697,  -105.25981970,  -105.45768990,  -105.69570612,  -105.81569229,  -105.37345137,  -105.32516163,
			 -105.18453534,  -104.77552447,  -105.66372894,  -104.74544305,  -105.20104130,  -105.62118650,  -105.94380708,  -105.47714940,
			 -105.19874793,  -105.68525828,  -105.25523303,  -106.00921372,  -105.42817906,  -105.37039116,  -104.92177477,  -105.12020005,
			 -105.08078506,  -105.10353405,  -104.54518056,  -105.42372788,  -104.36276160,  -105.57048219,  -104.25304132,  -105.44707834,
			 -105.85164970,  -105.31005854,  -105.97179071,  -104.91206281,  -105.69003208,  -105.51827466,  -104.98031915,  -105.57272140,
			 -104.69168806,  -104.40878076,  -104.86724947,  -104.93314076,  -105.42707029,  -105.34975959,  -104.90747149,  -105.62328258,
			 -105.75073814,  -104.86824433,  -105.33459151,  -104.57528847,  -105.69902440,  -104.64177276,  -105.30313879,  -105.00371520,
			 -105.54711433,  -104.08473685,  -105.57243620,  -104.97243129,  -104.69473109,  -105.43236988,  -105.00083400,  -105.71047073,
			 -104.82613769,  -105.57488855,  -104.05917939,  -105.31305372,  -105.41318547,  -104.32261347,  -104.59880840,  -105.91694962,
			 -105.33273335,  -105.42512567,  -104.14833099,  -105.78969791,  -104.52078892,  -104.37539318,  -106.21454791,  -105.17643874,
			 -104.41988197,  -104.65436707,  -106.16331996,  -105.40771477,  -104.77866487,  -104.83518526,  -104.09812772,  -105.19467591,
			 -104.81302279,  -105.06326489,  -104.85752660,  -104.47460566,  -105.48509406,  -103.97419265,  -104.72747351,  -106.26512238,
			 -105.30567872,  -104.62741587,  -104.39179561,  -106.19581383,  -104.90730321,  -104.58659249,  -104.00093258,  -105.29750442,
			 -105.16157435,  -104.44814917,  -104.90144203,  -104.89764416,  -106.10554281,  -105.94295558,  -105.92238245,  -103.77998285,
			 -105.24487515,  -105.04573505,  -104.12203689,  -104.32154083,  -104.94327460,  -104.89193137,  -104.46913816,  -104.62824592,
			 -105.26879510,  -105.54998474,  -105.53386738,  -105.73441590,  -105.25534260,  -105.67282977,  -105.41589873,  -103.86876868,
			 -105.21164829,  -105.41910383,  -104.38504775,  -105.31640073,  -106.02276866,  -105.95464693,  -104.61106233,  -104.99476755,
			 -106.19848688,  -104.55644170,  -105.90533052,  -105.24132852,  -105.93190072,  -105.13953344,  -105.09304051,  -103.82124776,
			 -105.03716091,  -104.52429798,  -104.86147576,  -104.53266535,  -103.67841444,  -104.67901055,  -104.68040140,  -104.21567143,
			 -104.38923651,  -105.61715733,  -104.90348587,  -105.84652956,  -104.82413091,  -105.19331684,  -105.75282474,  -105.24102379,
			 -104.91834144,  -106.01403595,  -105.23197119,  -105.37357883,  -106.14771254,  -103.86535116,  -103.91555954,  -104.49899593,
			 -106.07429149,  -105.94805621,  -104.30300892,  -104.59024552,  -104.83172202,  -105.61273732,  -105.51815959,  -105.04124813,
			 -104.27625896,  -104.12444693,  -105.29785190,  -103.72015330,  -104.89585332,  -104.67573400,  -104.78906648,  -104.71149972
		};
		
		final double Epsilon = 1e-3;
		
		int n = minimizedConfs.size();
		for (int i=0; i<n; i++) {
			
			double energy = minimizedConfs.get(i).getEnergy();
			
			if (i < expectedEnergies.length) {
				
				double absErr = energy - expectedEnergies[i];
				if (absErr > Epsilon) {
					
					System.out.println(String.format("\tWARNING: low precision energy: i:%-3d  exp:%12.8f  obs:%12.8f  absErr:%12.8f",
						i, expectedEnergies[i], energy, absErr
					));
					
				} else if (absErr < -Epsilon) {
				
					System.out.println(String.format("\timproved energy: i:%-3d  exp:%12.8f  obs:%12.8f  improvement:%12.8f",
						i, expectedEnergies[i], energy, -absErr
					));
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
