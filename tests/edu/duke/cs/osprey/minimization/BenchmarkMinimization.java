package edu.duke.cs.osprey.minimization;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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

public class BenchmarkMinimization extends TestBase {
	
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
		
		// settings
		final int numConfs = 40;
		final int[] numThreadsList = { 1, 2, 4, 8 };
		final boolean useGpu = true;
		
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
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, search.confSpace, search.shellResidues);
		
		// get the energy matrix
		File ematFile = new File("/tmp/emat.benchmarkMinimization.dat");
		if (ematFile.exists()) {
			System.out.println("\nReading energy matrix...");
			search.emat = (EnergyMatrix)ObjectIO.readObject(ematFile.getAbsolutePath(), true);
		}
		if (search.emat == null) {
			System.out.println("\nComputing energy matrix...");
			search.emat = new SimpleEnergyMatrixCalculator(ecalc).calcEnergyMatrix(); 
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
		}
		
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
		
		// what do we expect the energies to be?
		double[] expectedEnergies = {
			-107.01471465433335,
			-107.14427781940432,
			-106.79145713231975,
			-106.92139365967053,
			-106.28769308885211,
			-107.11801397703762,
			-106.678922061133,
			-106.41908247351522,
			-106.89600279606412,
			-106.74468003314176,
			-106.45689550906734,
			-106.95533592350961,
			-106.43786020587382,
			-106.52194038276753,
			-106.39164304147317,
			-106.72599111242266,
			-106.10055039136107,
			-106.73371625732581,
			-106.45590591321724,
			-105.95183095436968,
			-106.40589873963415,
			-106.16492362939529,
			-106.50474766966147,
			-105.73721952508399,
			-105.82709235989205,
			-106.01931524286667,
			-106.23489219140693,
			-106.38648146252073,
			-106.23062800511192,
			-106.01634808496966,
			-106.51221811051244,
			-106.13624568652024,
			-105.51519524725856,
			-105.58610640846234,
			-106.16496203741923,
			-105.99987206614541,
			-105.72752981857495,
			-105.98691660302913,
			-106.18174504161567,
			-105.79283093027898
		};
		
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
		checkEnergies(expectedEnergies, minimizedConfs);
		
		// benchmark parallel minimization
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		
		for (int numThreads : numThreadsList) {
			
			tasks.start(numThreads);
			
			System.out.println("\nBenchmarking " + numThreads + " task thread(s)...");
			Stopwatch taskStopwatch = new Stopwatch().start();
			minimizedConfs = new ConfMinimizer().minimize(confs, efuncs, search.confSpace, tasks);
			taskStopwatch.stop();
			System.out.println(String.format("Speedup: %.2fx", (float)baseStopwatch.getTimeNs()/taskStopwatch.getTimeNs()));
			checkEnergies(expectedEnergies, minimizedConfs);
			
			tasks.stopAndWait(10000);
		}
	}

	private static void checkEnergies(double[] expectedEnergies, List<EnergiedConf> minimizedConfs) {
		
		final double Epsilon = 1e-6;
		
		int n = minimizedConfs.size();
		for (int i=0; i<n; i++) {
			
			double energy = minimizedConfs.get(i).getEnergy();
			double absErr = energy - expectedEnergies[i];
			double relErr = absErr/Math.abs(expectedEnergies[i]);
			if (relErr > Epsilon) {
				System.out.println(String.format("\tWARNING: low precision energy: exp:%12.8f  obs:%12.8f  absErr:%12.8f  relErr:%12.8f",
					expectedEnergies[i], energy,
					absErr, relErr
				));
			}
		}
	}
}
