package edu.duke.cs.osprey.minimization;

import java.io.File;
import java.util.ArrayList;
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
import edu.duke.cs.osprey.energy.ForcefieldInteractionsGenerator;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
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
		final int numConfs = 256;//1024;//8192;
		
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
		
		benchmarkSerial(search, confs);
		benchmarkParallel(search, confs);
	}

	private static void benchmarkSerial(SearchProblem search, List<ScoredConf> confs)
	throws Exception {

		ForcefieldInteractionsGenerator ffintergen = new ForcefieldInteractionsGenerator();
		Factory<ForcefieldInteractions,Molecule> interactionsFactory = (mol) -> ffintergen.makeFullConf(search.confSpace, search.shellResidues, mol);
		ForcefieldParams ffparams = makeDefaultFFParams();
		
		System.out.println("\nbenchmarking CPU original...");
		Stopwatch cpuOriginalStopwatch = new Stopwatch();
		{
			Factory<Minimizer,MoleculeModifierAndScorer> minimizers = (mof) -> new CCDMinimizer(mof, false);
			Factory<EnergyFunction,Molecule> efuncs = (mol) -> EnvironmentVars.curEFcnGenerator.fullConfEnergy(search.confSpace, search.shellResidues, mol);
			cpuOriginalStopwatch.start();
			List<EnergiedConf> minimizedConfs = new ConfMinimizer(minimizers).minimize(confs, efuncs, search.confSpace);
			cpuOriginalStopwatch.stop();
			System.out.println(String.format("precise timing: %s, ops: %.1f", cpuOriginalStopwatch.getTime(TimeUnit.MILLISECONDS), confs.size()/cpuOriginalStopwatch.getTimeS()));
			checkEnergies(minimizedConfs);
		}
		
		System.out.println("\nbenchmarking CPU conf minimizer...");
		Stopwatch cpuSimpleStopwatch = benchmark(new CpuConfMinimizer(1, ffparams, interactionsFactory, search.confSpace), confs, cpuOriginalStopwatch);
		
		System.out.println("\nbenchmarking OpenCL simple...");
		benchmark(new GpuConfMinimizer(GpuConfMinimizer.Type.OpenCL, 1, 1, ffparams, interactionsFactory, search.confSpace), confs, cpuSimpleStopwatch);
		
		System.out.println("\nbenchmarking Cuda CCD...");
		benchmark(new GpuConfMinimizer(GpuConfMinimizer.Type.Cuda, 1, 1, ffparams, interactionsFactory, search.confSpace), confs, cpuSimpleStopwatch);
	}
	
	private static void benchmarkParallel(SearchProblem search, List<ScoredConf> confs)
	throws Exception {
		
		// settings
		final int[] numThreadsList = { 1 };//, 2, 4, 8 };
		final int[] numStreamsList = { 1, 2, 4, 8, 16 };//, 32, 64, 128, 256 };
		int maxNumStreams = numStreamsList[numStreamsList.length - 1];
		
		ForcefieldInteractionsGenerator ffintergen = new ForcefieldInteractionsGenerator();
		Factory<ForcefieldInteractions,Molecule> interactionsFactory = (mol) -> ffintergen.makeFullConf(search.confSpace, search.shellResidues, mol);
		ForcefieldParams ffparams = makeDefaultFFParams();
		
		List<EnergiedConf> minimizedConfs;
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		
		// benchmark cpu
		Stopwatch oneCpuStopwatch = null;
		for (int numThreads : numThreadsList) {
			
			tasks.start(numThreads);
			
			System.out.println("\nBenchmarking " + numThreads + " thread(s) with CPU efuncs...");
			Stopwatch stopwatch = benchmark(new CpuConfMinimizer(numThreads, ffparams, interactionsFactory, search.confSpace), confs, null);
			if (oneCpuStopwatch == null) {
				oneCpuStopwatch = stopwatch;
			}
			
			tasks.stopAndWait(10000);
		}
		
		// benchmark opencl
		Stopwatch oneOpenCLStopwatch = null;
		for (int numStreams : numStreamsList) {
			
			tasks.start(numStreams);
			
			System.out.println("\nBenchmarking " + numStreams + " stream(s) with OpenCL efuncs...");
			benchmark(new GpuConfMinimizer(GpuConfMinimizer.Type.OpenCL, 1, numStreams, ffparams, interactionsFactory, search.confSpace), confs, oneCpuStopwatch);
			
			tasks.stopAndWait(10000);
		}
		
		// benchmark cuda
		Stopwatch oneCudaStopwatch = null;
		for (int numStreams : numStreamsList) {
			
			tasks.start(numStreams);
			
			System.out.println("\nBenchmarking " + numStreams + " stream(s) with Cuda CCD minimizer...");
			benchmark(new GpuConfMinimizer(GpuConfMinimizer.Type.Cuda, 1, numStreams, ffparams, interactionsFactory, search.confSpace), confs, oneCpuStopwatch);
			
			tasks.stopAndWait(10000);
		}
	}
	
	private static Stopwatch benchmark(SpecializedConfMinimizer minimizer, List<ScoredConf> confs, Stopwatch referenceStopwatch)
	throws Exception {
		
		Stopwatch stopwatch = new Stopwatch().start();
		List<EnergiedConf> minimizedConfs = minimizer.minimize(confs);
		stopwatch.stop();
		
		System.out.print(String.format("precise timing: %s, ops: %.1f", stopwatch.getTime(TimeUnit.MILLISECONDS), confs.size()/stopwatch.getTimeS()));
		
		if (referenceStopwatch != null) {
			System.out.println(String.format(", speedup: %.2fx",
				(double)referenceStopwatch.getTimeNs()/stopwatch.getTimeNs()
			));
		} else {
			System.out.println();
		}
		
		minimizer.cleanup();
		
		checkEnergies(minimizedConfs);
		
		return stopwatch;
	}
	
	private static void checkEnergies(List<EnergiedConf> minimizedConfs) {
		
		// what do we expect the energies to be?
		final double[] expectedEnergies = {
		  -89.40965922,   -89.10803783,   -89.80947890,   -88.87688224,   -89.12763651,   -89.50404206,   -88.39800225,   -88.88968938,
		  -88.91538573,   -88.37400350,   -88.72521732,   -88.95813013,   -88.71957443,   -89.13542393,   -88.39342982,   -88.61512253,
		  -87.67170876,   -87.75762493,   -88.82437753,   -87.49111955,   -88.39267148,   -88.73093938,   -88.16624474,   -88.47911298,
		  -88.77040379,   -88.59812445,   -88.63729374,   -88.83029528,   -88.07022384,   -87.82115594,   -89.15034617,   -89.52776262,
		  -88.08922551,   -87.24538989,   -87.66296341,   -89.47261229,   -88.53548693,   -88.21416862,   -87.18239056,   -88.37126489,
		  -88.91533055,   -88.95432276,   -88.34024189,   -88.53617041,   -87.76065876,   -87.75246826,   -89.32887293,   -89.12214183,
		  -87.53435849,   -89.30674536,   -88.86121108,   -88.00498514,   -89.24745408,   -86.93536186,   -87.83485265,   -89.18378421,
		  -87.60530136,   -87.88059458,   -88.99239407,   -89.00570101,   -88.47514883,   -88.62549053,   -89.05482774,   -88.65430730,
		  -86.93662833,   -87.47811578,   -88.24904718,   -89.04046547,   -87.84106817,   -88.84132697,   -87.69399696,   -89.04363586,
		  -87.26109965,   -88.18510040,   -88.09667033,   -88.82233397,   -87.84429276,   -87.64489462,   -88.73508648,   -88.54431967,
		  -88.64763201,   -87.35367055,   -87.99012191,   -87.03408806,   -88.70440919,   -88.71415398,   -86.45949518,   -88.18726165,
		  -88.71381901,   -87.25162092,   -88.76650007,   -87.11761940,   -87.53086004,   -87.03019605,   -88.17094626,   -87.72276499,
		  -88.55346477,   -88.32758335,   -87.88605048,   -88.50483045,   -88.49351431,   -87.03813015,   -88.52753692,   -88.24436437,
		  -88.38377672,   -86.77330986,   -87.70085515,   -86.21378289,   -86.75358935,   -87.93605536,   -87.53739303,   -87.92238349,
		  -87.89511194,   -87.92618934,   -88.43244835,   -88.21239388,   -87.81246942,   -88.34964714,   -87.87499001,   -87.76636371,
		  -88.12503931,   -88.26555180,   -87.91955816,   -88.02610847,   -88.02430170,   -86.54615996,   -86.37892894,   -87.86968832,
		  -88.01315913,   -86.44454184,   -87.49567284,   -86.47629085,   -87.22640745,   -87.64994959,   -87.28919389,   -87.62895923,
		  -88.07505854,   -87.63801758,   -87.73628555,   -87.50320562,   -87.92081107,   -87.57626507,   -87.46741755,   -87.81516769,
		  -86.80913821,   -87.71664516,   -87.52480563,   -87.01066694,   -87.15025057,   -87.67471920,   -87.50873378,   -86.99428611,
		  -87.44875844,   -87.11867338,   -87.11426717,   -87.57375258,   -85.83146049,   -87.27720380,   -87.22563460,   -87.22713221,
		  -87.26064262,   -87.36576222,   -87.63173129,   -86.94828655,   -85.72236130,   -86.43930315,   -86.50554000,   -86.78037819,
		  -86.84000603,   -86.73110482,   -86.47304251,   -87.09353112,   -86.97871536,   -86.92161300,   -86.95181776,   -87.32759091,
		  -86.99452455,   -86.75702941,   -86.80736754,   -86.51632985,   -86.16466245,   -86.49014396,   -86.64316781,   -86.50261048,
		  -86.21235801,   -85.15617600,   -86.24234369,   -86.19326651,   -86.22862756,   -86.59775996,   -85.15255267,   -85.44153874,
		  -84.85099862,   -84.84617013,   -85.36785007,   -85.77339414,   -85.48292183,   -85.89434570,   -85.06227326,   -81.83221781,
		  -84.13343433,   -85.51132512,   -85.72983349,   -85.45257765,   -85.58801992,   -81.56680323,   -84.11896582,   -81.66978110,
		  -85.23320020,   -85.60141112,   -85.20586212,   -85.08154775,   -85.42391821,   -84.42927395,   -85.14686378,   -84.34419261,
		  -84.51351294,   -84.36367620,   -85.08003613,   -85.26103550,   -83.99720548,   -84.90237721,   -81.40550417,   -84.92652513,
		  -85.29660553,   -84.77621584,   -84.12384649,   -84.87167288,   -84.83449357,   -84.20914139,   -84.05879138,   -80.80884661,
		  -84.95585102,   -84.77254833,   -84.59822010,   -83.69101702,   -84.48839348,   -84.70707290,   -83.66418817,   -84.42932935,
		  -84.53021230,   -84.21010357,   -80.64689458,   -84.05838464,   -83.36037390,   -84.56809240,   -83.40692145,   -83.34103648,
		  -83.48031325,   -84.05700481,   -82.97432081,   -84.22780447,   -83.86858440,   -81.33965869,   -81.48681784,   -80.84873725,
		  -83.80115810,   -81.07571123,   -81.22193969,   -80.58561034,   -82.63094732,   -80.88584902,   -80.36449298,   -80.31331258,
		  -80.62236226,   -80.31674759,   -80.46338374,   -80.10152835,   -80.31088475,   -80.05014031,   -79.82585967,   -80.04816179,
		  -79.76886632,   -79.86249019,   -79.91703878,   -79.84850712,   -79.69191601,   -79.50500340,   -79.34116898,   -79.28969214,
		  -79.65392504,   -79.58594391,   -79.39012531,   -79.28785390,   -79.86225120,   -79.20646206,   -79.41059504,   -79.52178862,
		  -80.62886795,   -79.29490628,   -79.92881660,   -80.24364292,   -79.87516440,   -78.74583138,   -79.57900245,   -79.29103273,
		  -79.55894297,   -78.90200106,   -79.43596606,   -79.69202240,   -80.34000565,   -78.98620969,   -79.26067787,   -78.89430514,
		  -78.82564345,   -79.94575988,   -80.28542116,   -78.65686397,   -79.56563037,   -79.34394072,   -79.13596259,   -79.46531304,
		  -79.47974246,   -78.95655551,   -78.44429941,   -80.25679342,   -79.07421015,   -78.70918023,   -79.99750020,   -78.97231425,
		  -79.30563483,   -78.90158192,   -79.38257407,   -79.64896538,   -79.25333460,   -78.74649987,   -78.82573967,   -78.16695992,
		  -79.28758349,   -80.33566845,   -78.77885606,   -78.26346076,   -79.96885967,   -78.83870041,   -80.27903740,   -79.17799726,
		  -79.78964493,   -78.95427580,   -79.47190649,   -78.84325132,   -79.68387833,   -78.51164101,   -79.06070234,   -78.57068159,
		  -79.00080529,   -80.05370507,   -79.87824485,   -79.63988059,   -79.97319863,   -79.67617904,   -78.17395780,   -79.99026543,
		  -79.94395555,   -79.49661072,   -78.65660150,   -78.22222717,   -78.72800061,   -77.88774808,   -78.48429454,   -79.37035027,
		  -78.52757648,   -79.58027541,   -79.39407066,   -78.98354317,   -79.76236163,   -79.68129210,   -78.11651297,   -79.74491939,
		  -78.26166921,   -78.69628353,   -79.27750728,   -79.59715145,   -78.48850074,   -79.68477228,   -79.34974528,   -79.38801373,
		  -79.54338045,   -78.45202151,   -79.68101087,   -79.46210172,   -77.94776497,   -79.38708894,   -77.58390206,   -78.33430213,
		  -78.81787204,   -79.07964272,   -79.23988130,   -79.48269786,   -79.47438668,   -77.68065143,   -79.38980914,   -79.45661733,
		  -77.85386862,   -79.30011890,   -78.17835808,   -78.80633669,   -79.28917034,   -78.68024972,   -79.08468344,   -78.15020911,
		  -78.48344547,   -78.84540367,   -78.82454116,   -79.14570943,   -79.16760271,   -79.10264335,   -78.51343534,   -78.95228377,
		  -78.42731931,   -78.57830815,   -78.70478177,   -78.30587458,   -79.20136302,   -77.37631964,   -78.92776774,   -78.88720647,
		  -78.50643468,   -78.90604749,   -78.99416954,   -74.54888607,   -78.69746449,   -78.01955353,   -77.75519914,   -76.92053687,
		  -78.44321632,   -79.00131149,   -78.16342830,   -78.62535513,   -77.53971080,   -77.96176182,   -79.47916468,   -78.66524690,
		  -78.46423318,   -77.25586724,   -78.66657575,   -78.60848948,   -78.78329653,   -77.99492238,   -76.84905228,   -77.85280229,
		  -78.68986702,   -78.76531470,   -78.61801298,   -78.19719260,   -78.01604913,   -77.83373023,   -78.95259114,   -74.24465106,
		  -77.71622484,   -78.09813674,   -76.61692455,   -78.46637043,   -79.27612893,   -79.07169221,   -77.00920492,   -79.27499303,
		  -77.45670554,   -78.26085957,   -78.35582401,   -77.41476677,   -77.28349575,   -78.43392534,   -78.71791163,   -77.77842765,
		  -79.13194022,   -78.95576690,   -78.71887463,   -77.70760446,   -76.64231125,   -78.30961850,   -78.75765063,   -78.79581614,
		  -76.97945437,   -78.61025619,   -77.89151533,   -77.99854602,   -78.12499535,   -78.27761881,   -77.19496437,   -78.85710389,
		  -77.29353292,   -77.27420094,   -78.66416885,   -76.94535655,   -77.60526579,   -77.89644094,   -77.95006498,   -77.47451345,
		  -73.51512687,   -76.98031180,   -75.88155792,   -78.31160040,   -77.79163648,   -77.58599056,   -77.69000414,   -77.63384364,
		  -77.85825361,   -77.85828196,   -77.51595483,   -78.35614368,   -76.64013930,   -76.65050805,   -76.98477854,   -76.99379435,
		  -76.51892073,   -77.73037220,   -77.30188672,   -76.25029067,   -77.40252469,   -77.57689474,   -77.18615208,   -77.21123275,
		  -77.63937792,   -77.53070149,   -77.89059457,   -77.79128338,   -76.68777984,   -77.53951201,   -76.96669197,   -76.85781434,
		  -75.91489610,   -75.92281719,   -77.28944615,   -77.44027481,   -75.60854697,   -76.48216230,   -76.20897841,   -77.04211202,
		  -76.99317190,   -77.02549849,   -75.83901191,   -75.97062537,   -77.39785191,   -76.70313134,   -75.86926260,   -72.24506197,
		  -76.23851497,   -75.90345129,   -75.53457681,   -76.26858735,   -76.57363453,   -76.36439037,   -76.39687105,   -75.56430699,
		  -71.97985452,   -76.28326683,   -76.12795247,   -74.88148866,   -72.04422922,   -75.41439652,   -75.87472195,   -76.23319818,
		  -75.96315732,   -75.82337949,   -76.05849740,   -75.82220376,   -76.06468041,   -75.84886015,   -75.18686161,   -71.78018987,
		  -75.72111348,   -75.11062152,   -74.91111243,   -75.56813137,   -75.92840918,   -74.80621705,   -75.51811650,   -74.91918325,
		  -75.68064019,   -74.84686268,   -75.62192598,   -70.31024559,   -71.22206307,   -75.75950773,   -75.54140527,   -75.41698980,
		  -75.24602633,   -75.34181064,   -75.12125757,   -75.10487665,   -75.31768390,   -70.00501475,   -70.14355174,   -71.02157444,
		  -74.85178017,   -74.38139259,   -74.80037554,   -75.20005979,   -75.65008469,   -71.62934120,   -74.82600798,   -75.03162051,
		  -69.83331727,   -70.24688389,   -75.26713457,   -75.48515083,   -74.68749977,   -71.47749764,   -71.60707513,   -75.20716032,
		  -71.12508721,   -74.58876495,   -69.27169230,   -71.46820409,   -70.04164042,   -69.98350181,   -74.98756458,   -70.38501594,
		  -75.36873073,   -74.83654318,   -70.95896757,   -74.18487364,   -71.21380907,   -71.34234155,   -74.28071400,   -70.86198213,
		  -71.11742311,   -74.11929341,   -75.02876409,   -74.83466895,   -74.66929396,   -73.75127529,   -69.77877519,   -69.34410246,
		  -69.11254742,   -70.66955835,   -70.78433763,   -74.60209960,   -70.61122064,   -70.85391685,   -69.70508414,   -69.11477335,
		  -69.39988295,   -69.22321956,   -73.43171862,   -69.76264545,   -70.34823183,   -70.45476793,   -70.58410141,   -69.44179033,
		  -70.10233496,   -68.91300619,   -69.92336084,   -69.98001992,   -69.13723711,   -68.65509978,   -69.01880864,   -69.46640343,
		  -68.51250021,   -70.03059074,   -68.67832152,   -70.11708180,   -69.37050618,   -70.09417214,   -68.53034874,   -68.58836719,
		  -71.13773421,   -71.28357017,   -68.43250772,   -70.36235434,   -69.57816632,   -70.64730601,   -68.86256567,   -70.69410240,
		  -69.58801411,   -68.46272109,   -68.97947212,   -69.31033612,   -68.68258987,   -70.09834763,   -70.55385947,   -68.15854931,
		  -67.95675967,   -68.37724312,   -69.31376627,   -69.96186526,   -70.68339977,   -67.69901227,   -69.36438952,   -70.22777195,
		  -69.14317059,   -69.03156366,   -69.71398187,   -70.28347880,   -69.06703244,   -70.16240573,   -67.57530251,   -70.11012238,
		  -68.79326333,   -69.07725732,   -68.88691905,   -69.20231378,   -69.19314803,   -68.97000578,   -68.88186806,   -68.77823391,
		  -69.55146063,   -68.51690905,   -68.78008849,   -68.58519535,   -69.78267426,   -67.37297029,   -70.10887895,   -67.98934419,
		  -68.58878746,   -68.51749979,   -68.68261904,   -69.15008128,   -68.29886801,   -68.36896279,   -68.28723470,   -69.04772809,
		  -67.69423043,   -68.38458738,   -69.56593665,   -68.88327670,   -68.77390703,   -68.44845853,   -68.75062662,   -68.21900701,
		  -69.06394455,   -68.17706323,   -68.96890429,   -68.08734748,   -68.65431258,   -68.87630270,   -68.96468773,   -67.35089158,
		  -67.71433100,   -69.71534198,   -69.64701354,   -67.46336280,   -69.46104033,   -68.60548192,   -67.99052400,   -68.01573489,
		  -67.85148827,   -68.53294453,   -66.62593808,   -68.39641382,   -68.22379788,   -68.61368161,   -67.35527039,   -67.54984824,
		  -67.70023499,   -69.37651423,   -66.73390476,   -69.52783822,   -70.05426079,   -68.22167158,   -67.70779668,   -69.93353589,
		  -69.63563212,   -68.97658271,   -70.43030110,   -66.96823431,   -68.28566922,   -67.48235544,   -67.60554195,   -69.06010411,
		  -68.15361101,   -68.19906359,   -67.39245858,   -68.56745030,   -69.59730357,   -67.50342153,   -67.47498177,   -69.63912449,
		  -70.21574958,   -66.47497573,   -67.29732070,   -67.41244790,   -69.02707141,   -66.24567552,   -68.98792335,   -67.74436235,
		  -66.16899373,   -69.47708559,   -66.97380440,   -69.34797635,   -67.45580726,   -66.16813339,   -66.98370919,   -68.93294425,
		  -70.05966437,   -70.28776802,   -66.33177018,   -66.66830340,   -68.76831395,   -69.15388861,   -65.53823050,   -66.57381041,
		  -70.08302552,   -69.58827785,   -68.67278099,   -69.82907377,   -69.59451184,   -69.91170257,   -69.48693072,   -66.06287960,
		  -69.32884185,   -65.45256290,   -66.50934161,   -69.77673931,   -69.63071028,   -68.33655918,   -69.47980782,   -68.76829416,
		  -67.65953028,   -69.40041722,   -69.56619653,   -69.34385307,   -69.54762648,   -65.98431793,   -67.06828692,   -69.43351815,
		  -68.25310261,   -68.21863167,   -69.04353639,   -65.68206629,   -68.58445857,   -66.68528278,   -66.90338651,   -67.45275132,
		  -66.62568115,   -65.67987268,   -69.09229384,   -65.14413219,   -65.52728198,   -66.40198612,   -65.37777745,   -66.25494169,
		  -68.07084046,   -65.60455121,   -57.57701299,   -68.69382121,   -68.27129730,   -64.31547631,   -67.78535564,   -65.53795824,
		  -66.68429122,   -64.67447542,   -66.53916211,   -64.84078673,   -66.22609873,   -65.15685287,   -65.22255184,   -68.43088836,
		  -64.67635154,   -65.45100039,   -64.96116167,   -66.19916155,   -67.78195075,   -65.83985515,   -64.37130229,   -57.44337377,
		  -67.05088849,   -64.65899960,   -60.45474198,   -68.19990133,   -64.37273332,   -65.77255926,   -64.60434092,   -64.50440648,
		  -64.11037846,   -67.65799284,   -60.19057198,   -67.76384543,   -66.70070484,   -63.64080994,   -63.64288611,   -67.28269512,
		  -61.44463689,   -66.74776037,   -61.19085814,   -59.43196587,   -61.18244184,   -57.09970967,   -61.08905015,   -57.21877495,
		  -60.92882214,   -65.68199736,   -56.61528813,   -60.82656781,   -56.64548521,   -65.96498436,   -65.60695429,   -60.42128200,
		  -56.13775754,   -66.45934391,   -56.06771705,   -62.04304763,   -65.62530908,   -60.16759911,   -66.02470434,   -60.06599032,
		  -66.12017248,   -56.06998090,   -65.88295702,   -61.84329191,   -60.66998691,   -65.18243195,   -65.62953250,   -66.00098006,
		  -65.57881250,   -55.51909702,   -65.83286743,   -65.60396512,   -65.48850920,   -55.67644571,   -61.83548958,   -55.60347038,
		  -60.83970890,   -60.17406201,   -65.38999664,   -60.08133603,   -60.23629925,   -60.78436216,   -60.26524862,   -59.90782666,
		  -61.49715618,   -60.04354614,   -61.27615692,   -61.40483664,   -61.50047219,   -60.92404321,   -61.00455773,   -59.83998677,
		  -60.86694892,   -60.76098078,   -61.19758730,   -60.91146740,   -59.90958930,   -60.91529750,   -60.00426574,   -59.32286699,
		  -59.51473847,   -60.98895692,   -58.89244958,   -60.97374510,   -59.94554184,   -60.40949425,   -59.61438822,   -60.47430560,
		  -59.50376503,   -59.45800255,   -59.42983367,   -59.19881812,   -59.82064508,   -58.14041634,   -60.51855514,   -58.68536428,
		  -58.77369239,   -58.38119802,   -58.78635461,   -58.51631092,   -59.85247428,   -58.48074401,   -59.21693818,   -59.53861918,
		  -58.39164951,   -59.02328569,   -57.92505803,   -59.86657491,   -59.44447681,   -55.48489376,   -59.26132945,   -59.41622466,
		  -59.24794429,   -59.51293981,   -57.98292064,   -59.60379479,   -58.34995519,   -58.95714741,   -58.21728925,   -59.28352261,
		  -58.95679854,   -57.96441976,   -59.11328631,   -59.03848625,   -58.22331841,   -57.64227054,   -58.96012219,   -59.37283456,
		  -59.16734742,   -58.03348935,   -58.99061737,   -58.87029357,   -57.70209323,   -58.12167591,   -58.93688159,   -58.81782080,
		  -58.77334165,   -58.65151113,   -58.35005738,   -58.13380927,   -58.45294877,   -58.28242302,   -56.82520146,   -57.78098115,
		  -57.38541228,   -57.46575604,   -57.87693660,   -57.04498601,   -57.37113502,   -56.77887326,   -57.44247949,   -57.53792515,
		  -56.22929481,   -57.30097985,   -57.04684599,   -56.99667350,   -47.99263145,   -56.35468674,   -56.99492052,   -57.17097056,
		  -57.00271155,   -47.81064781,   -56.65859473,   -59.44406411,   -51.29071991,   -56.55997459,   -59.13920086,   -55.73970607
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
