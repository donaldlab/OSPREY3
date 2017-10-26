package edu.duke.cs.osprey.minimization;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.FFInterGen;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ResPairCache;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.cuda.Gpus;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.TimeFormatter;
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
		resFlex.sortPositions();
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
			"test", "examples/1CC8/1CC8.ss.pdb", 
			resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null,
			false, new ArrayList<>()
		);
		
		// calc the energy matrix
		File ematFile = new File("/tmp/benchmarkMinimization.emat.dat");
		if (ematFile.exists()) {
			search.emat = (EnergyMatrix)ObjectIO.readObject(ematFile.getAbsolutePath(), false);
		} else {
			search.emat = new SimpleEnergyMatrixCalculator.Cpu(2, makeDefaultFFParams(), search.confSpace, search.shellResidues).calcEnergyMatrix();
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
		}
		
		// also make a simple conf space
		Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
		strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
		strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
		strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
		strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
		strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
		strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
		strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
		SimpleConfSpace simpleConfSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		assertConfSpacesMatch(search.confSpace, simpleConfSpace);
		
		// settings
		final int numConfs = 1024*2;//32;//256;//512;//1024;
		
		// get a few arbitrary conformations
		System.out.println("getting confs...");
		search.pruneMat = new PruningMatrix(search.confSpace, 1000);
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setMPLP(new ConfAStarTree.MPLPBuilder()
				.setNumIterations(4)
			).build();
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
		
		// get all the confs
		List<ScoredConf> allConfs = new ArrayList<>(confs);
		while (true) {
			ScoredConf conf = tree.nextConf();
			if (conf == null) {
				break;
			}
			allConfs.add(conf);
		}
		
		System.out.println("benchmarking...");
		
		//benchmarkSerial(search, simpleConfSpace, confs);
		//benchmarkParallel(search, simpleConfSpace, confs);
		//compareOneConf(search, confs);
		
		benchmarkOps(search, simpleConfSpace, allConfs);
	}

	private static void benchmarkSerial(SearchProblem search, SimpleConfSpace simpleConfSpace, List<ScoredConf> confs)
	throws Exception {

		Factory<ForcefieldInteractions,Molecule> interactionsFactory = (mol) -> FFInterGen.makeFullConf(search.confSpace, search.shellResidues, mol);
		ForcefieldParams ffparams = makeDefaultFFParams();
		
		/*
		System.out.println("\nbenchmarking CPU original...");
		Stopwatch cpuOriginalStopwatch = new Stopwatch();
		{
			Factory<Minimizer,MoleculeModifierAndScorer> minimizers = (mof) -> new CCDMinimizer(mof, false);
			Factory<EnergyFunction,Molecule> efuncs = (mol) -> EnvironmentVars.curEFcnGenerator.fullConfEnergy(search.confSpace, search.shellResidues, mol);
			cpuOriginalStopwatch.start();
			List<EnergiedConf> minimizedConfs = new CpuConfMinimizer.Builder(ffparams, interactionsFactory, search.confSpace).setMinimizers(minimizers).build().minimize(confs);
			cpuOriginalStopwatch.stop();
			System.out.println(String.format("precise timing: %s, ops: %.1f", cpuOriginalStopwatch.getTime(TimeUnit.MILLISECONDS), confs.size()/cpuOriginalStopwatch.getTimeS()));
			checkEnergies(minimizedConfs);
		}
		*/
		
		System.out.println("\nbenchmarking CPU conf minimizer...");
		//Stopwatch cpuSimpleStopwatch = benchmark(new CpuConfMinimizer.Builder(ffparams, interactionsFactory, search.confSpace).build(), confs, cpuOriginalStopwatch);
		Stopwatch cpuSimpleStopwatch = benchmark(new CpuConfMinimizer.Builder(ffparams, interactionsFactory, search.confSpace).build(), confs, null);
		
		/*
		System.out.println("\nbenchmarking Residue CPU...");
		fragecalc = new MinimizingFragmentEnergyCalculator.Builder(simpleConfSpace, ffparams)
			.setType(MinimizingFragmentEnergyCalculator.Type.Cpu)
			.setParallelism(Parallelism.makeCpu(1))
			.build();
		benchmark(new MinimizingConfEnergyCalculator.Builder(fragecalc).build(), confs, cpuSimpleStopwatch);
		
		System.out.println("\nbenchmarking OpenCL simple...");
		benchmark(new GpuConfMinimizer.Builder(ffparams, interactionsFactory, search.confSpace).setGpuInfo(GpuConfMinimizer.Type.OpenCL, 1, 1).build(), confs, cpuSimpleStopwatch);
		
		System.out.println("\nbenchmarking Cuda...");
		benchmark(new GpuConfMinimizer.Builder(ffparams, interactionsFactory, search.confSpace).setGpuInfo(GpuConfMinimizer.Type.Cuda, 1, 1).build(), confs, cpuSimpleStopwatch);
		
		System.out.println("\nbenchmarking Cuda CCD...");
		benchmark(new GpuConfMinimizer.Builder(ffparams, interactionsFactory, search.confSpace).setGpuInfo(GpuConfMinimizer.Type.CudaCCD, 1, 1).build(), confs, cpuSimpleStopwatch);
		
		System.out.println("\nbenchmarking Residue Cuda...");
		fragecalc = new MinimizingFragmentEnergyCalculator.Builder(simpleConfSpace, ffparams)
			.setType(MinimizingFragmentEnergyCalculator.Type.ResidueCuda)
			.setParallelism(Parallelism.makeGpu(1, 1))
			.build();
		benchmark(new MinimizingConfEnergyCalculator.Builder(fragecalc).build(), confs, cpuSimpleStopwatch);
		*/
		
		System.out.println("\nbenchmarking Residue Cuda CCD...");
		new EnergyCalculator.Builder(simpleConfSpace, ffparams)
			.setType(EnergyCalculator.Type.ResidueCudaCCD)
			.setParallelism(Parallelism.make(4, 1, 1))
			.use((ecalc) -> {
				benchmark(simpleConfSpace, ecalc, confs, cpuSimpleStopwatch);
			}
		);
	}
	
	private static void benchmarkParallel(SearchProblem search, SimpleConfSpace simpleConfSpace, List<ScoredConf> confs)
	throws Exception {
		
		// settings
		final int[] numThreadsList = { 1 };//, 2, 4, 8 };
		final int numGpus = 1;
		final int[] numStreamsPerGpuList = { 1, 2, 4, 8, 16, 32 };//, 64, 128 };//, 256 };
		
		Factory<ForcefieldInteractions,Molecule> interactionsFactory = (mol) -> FFInterGen.makeFullConf(search.confSpace, search.shellResidues, mol);
		ForcefieldParams ffparams = makeDefaultFFParams();
		
		// benchmark cpu
		Stopwatch stopwatch = null;
		/*
		for (int numThreads : numThreadsList) {
			
			System.out.println(String.format("\nBenchmarking %2d thread(s) with %30s...", numThreads, "CPU efuncs"));
			ConfMinimizer minimizer = new CpuConfMinimizer.Builder(ffparams, interactionsFactory, search.confSpace)
				.setNumThreads(numThreads)
				.build();
			if (stopwatch == null) {
				stopwatch = benchmark(minimizer, confs, null);
			} else {
				benchmark(minimizer, confs, null);
			}
		}
		*/
		Stopwatch oneCpuStopwatch = stopwatch;
		
		/* benchmark cpu residue efuncs
		for (int numThreads : numThreadsList) {
			
			System.out.println(String.format("\nBenchmarking %3d thread(s) with %30s...", numThreads, "CPU residue efuncs"));
			new MinimizingFragmentEnergyCalculator.Builder(simpleConfSpace, ffparams)
				.setType(MinimizingFragmentEnergyCalculator.Type.Cpu)
				.setParallelism(Parallelism.makeCpu(numThreads))
				.use((fragEcalc) -> {
					MinimizingConfEnergyCalculator confEcalc = new MinimizingConfEnergyCalculator.Builder(fragEcalc).build();
					benchmark(confEcalc, confs, oneCpuStopwatch);
				});
		}
		*/
		
		/*
		// benchmark opencl
		for (int numStreams : numStreamsList) {
			System.out.println(String.format("\nBenchmarking %3d stream(s) with %30s...", numStreams, "OpenCL efuncs"));
			ConfMinimizer minimizer = new GpuConfMinimizer.Builder(ffparams, interactionsFactory, search.confSpace)
				.setGpuInfo(GpuConfMinimizer.Type.OpenCL, 1, numStreams)
				.build();
			benchmark(minimizer, confs, oneCpuStopwatch);
		}
		*/
		
		// benchmark cuda
		for (int numStreams : numStreamsPerGpuList) {
			System.out.println(String.format("\nBenchmarking %3d stream(s) with %30s...", numStreams, "Cuda efuncs"));
			ConfMinimizer minimizer = new GpuConfMinimizer.Builder(ffparams, interactionsFactory, search.confSpace)
				.setGpuInfo(GpuConfMinimizer.Type.Cuda, numGpus, numStreams)
				.build();
			benchmark(minimizer, confs, oneCpuStopwatch);
		}
		
		// benchmark cuda ccd
		for (int numStreams : numStreamsPerGpuList) {
			System.out.println(String.format("\nBenchmarking %3d stream(s) with %30s...", numStreams, "Cuda CCD minimizer"));
			ConfMinimizer minimizer = new GpuConfMinimizer.Builder(ffparams, interactionsFactory, search.confSpace)
				.setGpuInfo(GpuConfMinimizer.Type.CudaCCD, numGpus, numStreams)
				.build();
			benchmark(minimizer, confs, oneCpuStopwatch);
		}
		
		// benchmark residue cuda
		for (int numStreams : numStreamsPerGpuList) {
			System.out.println(String.format("\nBenchmarking %3d stream(s) with %30s...", numStreams, "Cuda residue efuncs"));
			new EnergyCalculator.Builder(simpleConfSpace, ffparams)
				.setType(EnergyCalculator.Type.ResidueCuda)
				.setParallelism(Parallelism.make(4, numGpus, numStreams))
				.use((ecalc) -> {
					benchmark(simpleConfSpace, ecalc, confs, oneCpuStopwatch);
				});
		}
		
		// benchmark residue cuda ccd
		for (int numStreams : numStreamsPerGpuList) {
			System.out.println(String.format("\nBenchmarking %3d stream(s) with %30s...", numStreams, "Cuda residue CCD minimizer"));
			new EnergyCalculator.Builder(simpleConfSpace, ffparams)
				.setType(EnergyCalculator.Type.ResidueCudaCCD)
				.setParallelism(Parallelism.make(4, numGpus, numStreams))
				.use((ecalc) -> {
					benchmark(simpleConfSpace, ecalc, confs, oneCpuStopwatch);
				});
		}
	}
	
	private static Stopwatch benchmark(ConfMinimizer minimizer, List<ScoredConf> confs, Stopwatch referenceStopwatch)
	throws Exception {
		
		minimizer.setReportProgress(true);
		
		Stopwatch stopwatch = new Stopwatch().start();
		List<EnergiedConf> minimizedConfs = minimizer.minimize(confs);
		stopwatch.stop();
		
		System.out.print(String.format("precise timing: %9s, ops: %5.2f", stopwatch.getTime(TimeUnit.MILLISECONDS), confs.size()/stopwatch.getTimeS()));
		
		if (referenceStopwatch != null) {
			System.out.println(String.format(", speedup: %6.2fx",
				(double)referenceStopwatch.getTimeNs()/stopwatch.getTimeNs()
			));
		} else {
			System.out.println();
		}
		
		minimizer.cleanup();
		
		checkEnergies(minimizedConfs);
		
		return stopwatch;
	}
	
	private static Stopwatch benchmark(SimpleConfSpace confSpace, EnergyCalculator ecalc, List<ScoredConf> confs, Stopwatch referenceStopwatch)
	throws Exception {
		
		ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();
		
		Stopwatch stopwatch = new Stopwatch().start();
		List<EnergiedConf> minimizedConfs = confEcalc.calcAllEnergies(confs, true);
		stopwatch.stop();
		
		System.out.print(String.format("precise timing: %9s, ops: %6.2f", stopwatch.getTime(TimeUnit.MILLISECONDS), confs.size()/stopwatch.getTimeS()));
		
		if (referenceStopwatch != null) {
			System.out.println(String.format(", speedup: %6.2fx",
				(double)referenceStopwatch.getTimeNs()/stopwatch.getTimeNs()
			));
		} else {
			System.out.println();
		}
		
		checkEnergies(minimizedConfs);
		
		return stopwatch;
	}
	
	private static void compareOneConf(SearchProblem search, List<ScoredConf> confs)
	throws Exception {
		
		ScoredConf conf = confs.get(125);
		
		Factory<ForcefieldInteractions,Molecule> ffinteractions = (mol) -> FFInterGen.makeFullConf(search.confSpace, search.shellResidues, mol);
		ForcefieldParams ffparams = makeDefaultFFParams();
		
		double originalEnergy = new CpuConfMinimizer.Builder(ffparams, ffinteractions, search.confSpace)
			.setMinimizers((mof) -> new CCDMinimizer(mof, false))
			.build().minimize(conf).getEnergy();
		
		double cpuEnergy = new CpuConfMinimizer.Builder(ffparams, ffinteractions, search.confSpace)
			.setMinimizers((mof) -> new SimpleCCDMinimizer(mof))
			.build().minimize(conf).getEnergy();
		
		double openclEnergy = new GpuConfMinimizer.Builder(ffparams, ffinteractions, search.confSpace)
			.setGpuInfo(GpuConfMinimizer.Type.OpenCL, 1, 1)
			.build().minimize(conf).getEnergy();
		
		double cudaEnergy = new GpuConfMinimizer.Builder(ffparams, ffinteractions, search.confSpace)
			.setGpuInfo(GpuConfMinimizer.Type.Cuda, 1, 1)
			.build().minimize(conf).getEnergy();
		
		double cudaCcdEnergy = new GpuConfMinimizer.Builder(ffparams, ffinteractions, search.confSpace)
			.setGpuInfo(GpuConfMinimizer.Type.CudaCCD, 1, 1)
			.build().minimize(conf).getEnergy();
		
		System.out.println(String.format("%21.16f", originalEnergy));
		System.out.println(String.format("%21.16f   %e", cpuEnergy, cpuEnergy - originalEnergy));
		System.out.println(String.format("%21.16f   %e", openclEnergy, openclEnergy - originalEnergy));
		System.out.println(String.format("%21.16f   %e", cudaEnergy, cudaEnergy - originalEnergy));
		System.out.println(String.format("%21.16f   %e", cudaCcdEnergy, cudaCcdEnergy - originalEnergy));
	}
	
	private static void checkEnergies(List<EnergiedConf> minimizedConfs) {
		
		// what do we expect the energies to be?
		final double[] expectedEnergies = {
			    -89.40966969380693,     -89.10792031498546,     -89.80959784194744,     -88.63999143595730,     -89.12813398453570,     -89.50404412354310,     -88.39619842054671,     -88.88944810225365,
			    -88.91539256575398,     -88.37401748233218,     -88.72521745746171,     -88.95852827530793,     -88.56492542985217,     -89.13542390896983,     -88.39342805726660,     -88.61512935924286,
			    -87.90968414029327,     -87.89000423010312,     -88.82437706142682,     -87.49112060974092,     -88.39275265803960,     -88.73114551918172,     -88.16537259873690,     -88.48123557923410,
			    -88.77040362062795,     -88.59813059337920,     -88.63729060259028,     -88.83029518512434,     -88.07022325053538,     -87.96393403364500,     -89.15033913943141,     -89.52692924559630,
			    -88.08922374497158,     -87.24538772522436,     -87.66296306373582,     -89.47261233439642,     -88.53548636087920,     -88.21415121135150,     -87.18239160797462,     -88.37126442878240,
			    -88.91533056304328,     -88.95432302575750,     -88.34009371038420,     -88.53616773986336,     -87.76065810577316,     -87.75248131792512,     -89.33135478542602,     -89.12228145797883,
			    -87.53441661042723,     -89.32471565270455,     -88.86123870156534,     -88.00546058189720,     -89.24541472451040,     -86.93535974296425,     -87.98785179871844,     -89.18378404733703,
			    -87.75560199933152,     -87.88060140510152,     -89.00171254523137,     -89.00570252622346,     -88.47514885622186,     -88.62549054867130,     -89.05482774046607,     -88.65432377503429,
			    -86.93663044950159,     -87.47811586993873,     -88.24899664839260,     -89.04419754206710,     -87.84105951054588,     -88.84132967727714,     -87.69399893847299,     -89.04235969695354,
			    -87.40533830584779,     -88.18513151332593,     -88.09667034751962,     -88.80872521464953,     -87.84429223080038,     -88.73509155911078,     -87.64489424433181,     -88.54657358501376,
			    -88.64814162469027,     -87.35366871460823,     -87.99012671191343,     -86.86384986410590,     -88.71371925646392,     -88.71789454575847,     -86.45949622431180,     -88.18726743622189,
			    -88.71134172319888,     -87.25162087938538,     -88.76650007186149,     -87.27974604769545,     -87.53085146547070,     -87.19438867495086,     -88.17093925109661,     -87.72276497741976,
			    -88.55346431246410,     -88.32681020819021,     -87.88608163341416,     -88.50405466420116,     -88.49351453892960,     -87.03812944503379,     -88.53063387303504,     -88.24438184069945,
			    -88.38437772566940,     -86.80662077101380,     -87.70085929866920,     -86.21378083766024,     -86.75366475185666,     -87.93605536725366,     -87.53739301960660,     -87.92240008983532,
			    -87.89511503886500,     -87.92736033689788,     -88.43244836981708,     -88.21984788323206,     -87.81246940565500,     -88.35298794228237,     -87.87499005959380,     -87.76636380712874,
			    -88.12500049661190,     -88.26555133854549,     -87.91957482925183,     -88.02602829918796,     -88.02430154732673,     -86.57899802649027,     -86.52638276331666,     -87.86983308364037,
			    -88.01375792432842,     -86.44461115537132,     -87.49567876076678,     -86.30136897748173,     -87.22640743882549,     -87.63192461850136,     -87.28919387771826,     -87.62863067216593,
			    -88.07505854152441,     -87.63919998467404,     -87.73602270463720,     -87.50320557241616,     -87.92081244613335,     -87.57626510288816,     -87.46741757103308,     -87.81475764520155,
			    -86.80912945601120,     -87.71656469237308,     -87.52481310839153,     -87.01067135174435,     -87.15028185313524,     -87.67464252156707,     -87.50875193830277,     -86.99428609082398,
			    -87.44849554308277,     -87.11908564421010,     -87.11426717413980,     -87.57375210674628,     -85.86045824866588,     -87.27719784540640,     -87.22563957372088,     -87.22713303464256,
			    -87.26064225630965,     -87.36568801101708,     -87.63168808674718,     -86.94826102252034,     -85.72243565594829,     -86.44304351031661,     -86.50553997605290,     -86.78037816675399,
			    -86.84000605366171,     -86.73110481810454,     -86.47304263273980,     -87.09353419769880,     -86.97870937794013,     -86.92161443514735,     -86.95181726029462,     -87.32758491152750,
			    -86.99444432483340,     -86.75676646880312,     -86.80736746223370,     -86.51633202193120,     -86.16466256668261,     -86.49015045866462,     -86.64308956159195,     -86.50261040703359,
			    -86.21236012169021,     -85.15617698487188,     -86.24233770927749,     -86.19326796635393,     -86.22862707702006,     -86.59771666514777,     -85.15255258570564,     -85.44153885361040,
			    -84.85099965287893,     -84.84617005072491,     -85.36784949502446,     -85.77339406014504,     -85.48292396860278,     -85.89433932360810,     -85.06227265031667,     -81.82767834392853,
			    -84.13343536504894,     -85.51133347249042,     -85.72984950163173,     -85.45257759792123,     -85.58801413359295,     -81.56525191833501,     -84.11896573075445,     -85.23317291647002,
			    -81.66970166044752,     -85.60142915591210,     -85.20587045884616,     -85.08154702198355,     -85.42393405357569,     -84.42934206695311,     -85.14686372754704,     -84.34419301336170,
			    -84.51350707740228,     -84.36341472002285,     -85.08003610901808,     -85.26102866183645,     -83.99720541655675,     -84.90237706407530,     -81.40535788267216,     -84.92649802710301,
			    -85.29662296604668,     -84.77621514198015,     -84.12366605501653,     -84.87166688040050,     -84.83449275823259,     -84.20913549835058,     -84.05852988831010,     -80.80760232168983,
			    -84.95584417230978,     -84.77254533733368,     -84.59822001681644,     -83.69101695360051,     -84.48838928302662,     -84.70708884373553,     -83.66418807865648,     -84.42932928711409,
			    -84.53021145472121,     -84.21007622678941,     -84.05838393044753,     -80.64676292174266,     -83.36037322071327,     -84.56811036929686,     -83.40696733404340,     -83.34077497171164,
			    -83.48030738885691,     -84.05700478476739,     -82.97432074404941,     -84.22779762114213,     -83.86858429026721,     -81.33965980291978,     -81.48681781383640,     -80.84873755603023,
			    -83.80115730411302,     -81.07571224229440,     -81.22193973918914,     -80.58561064645116,     -82.63100416096798,     -80.88584908052367,     -80.36521103582047,     -80.31331270948279,
			    -80.62236231457074,     -80.31674763597450,     -80.46338378496934,     -80.10224648296915,     -80.31088453128808,     -80.05014130117075,     -79.82585996904007,     -80.04816185160544,
			    -79.76886640628965,     -79.86249023513600,     -79.91703918327276,     -79.84966807825514,     -79.69191582855380,     -79.50500343035970,     -79.34188713013660,     -79.28969225231351,
			    -79.65392630027743,     -79.58710180012723,     -79.39012668969006,     -79.28785284214106,     -79.86225100800338,     -79.20646342683138,     -79.41069841369570,     -79.52172036969213,
			    -80.62886777918980,     -79.29490622521243,     -79.92933363542224,     -80.24364275372037,     -79.87516420138822,     -78.74583061254756,     -79.58205583130307,     -79.16917844788065,
			    -79.55894433368167,     -78.90200090332169,     -79.28563322361211,     -79.69105067204084,     -80.34000686953631,     -78.98620948135118,     -79.26067776031731,     -78.89430553616744,
			    -78.82680466271812,     -79.94575971145525,     -80.28542104864408,     -78.65686377613048,     -79.56563016416578,     -79.49504615535838,     -79.13618339797580,     -79.46808405844979,
			    -79.47974240405334,     -78.95656603817457,     -78.44450228448981,     -80.25679324289344,     -79.07299157052380,     -78.56388559345511,     -79.99750150016960,     -78.96770359205838,
			    -79.30534945484176,     -78.89926861765238,     -79.23917010232734,     -79.64896519618617,     -79.25333444160732,     -78.74979241750212,     -78.82573947520557,     -78.16695975480737,
			    -79.28758431603657,     -80.33431463717051,     -78.62020969616194,     -78.26346061509165,     -79.96885948804866,     -78.83764716166566,     -80.27903715356169,     -79.17801294079462,
			    -79.78964472879197,     -78.81163464113163,     -79.47190771840133,     -78.84325101217360,     -79.68387813340408,     -78.66378221420078,     -79.06070248120260,     -78.57040319900172,
			    -79.00080352150636,     -80.05356850460677,     -79.87808690017438,     -79.63988310173923,     -79.97269577889081,     -79.67610579987874,     -78.17393966854750,     -79.99026522187091,
			    -79.96158137971837,     -79.49661025862459,     -78.65623297340943,     -78.22222698654875,     -78.57682360046279,     -77.88774941030074,     -78.48430790726930,     -79.37035078579486,
			    -78.52757791691010,     -79.58027734969342,     -79.39407046129372,     -78.98077347105583,     -79.76236142295397,     -79.68129378243336,     -78.11965397370446,     -79.74491918070501,
			    -78.26140154558504,     -78.69628332189168,     -79.27750709753204,     -79.59715261754388,     -78.48849935421374,     -79.68430501024140,     -79.34974672482380,     -79.38793064415432,
			    -79.54338023314564,     -78.45202131830588,     -79.67960095785584,     -79.44880258520682,     -77.78761775697997,     -79.38709079750850,     -77.58390190775079,     -78.18832930594770,
			    -78.81787177572757,     -79.07964396145799,     -79.23988112328150,     -79.48035679162058,     -79.47438639394835,     -77.68065126043032,     -79.38981088546221,     -79.45661712115651,
			    -78.00628491531933,     -79.30011866325945,     -78.17835683679043,     -78.80633646849033,     -79.28917021177118,     -78.68025649259374,     -79.08468323373434,     -78.15020890219935,
			    -78.34786877632420,     -78.85057041471937,     -78.82495598859845,     -79.14495917670095,     -79.16760250476783,     -79.12261637219858,     -78.51343513186447,     -78.95228359933716,
			    -78.57691598197863,     -78.58570276903096,     -78.70478156377355,     -78.30587435393510,     -79.20136283670844,     -77.37631946133449,     -78.92776755607545,     -78.89514229245740,
			    -78.50643444016444,     -78.90604742941628,     -78.99367105452592,     -74.54888586874586,     -78.69739654253158,     -78.01964704590081,     -77.75519870875658,     -76.92053893535173,
			    -78.44317477457025,     -79.00131129569125,     -78.16342659384846,     -77.53944046959501,     -78.62535459098306,     -77.96176161055630,     -79.47918403292881,     -78.66524815472079,
			    -78.46423323379888,     -77.25586825897669,     -78.66657710549518,     -78.60851376536979,     -78.78329632132326,     -77.99492214957931,     -76.84905209674280,     -77.67697732054047,
			    -78.68986681617778,     -78.76531447722647,     -78.61801272909908,     -78.19719236022775,     -78.01604892925258,     -77.83372954280064,     -78.95258406070728,     -74.24465075615356,
			    -77.71622463745811,     -78.08149687202966,     -76.61678296626740,     -78.46637082277390,     -79.27612875006169,     -79.07242419332394,     -77.00920276466759,     -79.27445210710829,
			    -77.45670553939665,     -78.26085939089904,     -78.35581802456778,     -77.41476656463253,     -77.28352434809780,     -78.43392510627532,     -78.71791160402776,     -77.77842744398185,
			    -79.13579242513511,     -78.96347572616169,     -78.71887453245755,     -77.70760426162806,     -76.64231107576275,     -78.30961830609930,     -78.76753878510827,     -78.79637438985522,
			    -76.97948184127203,     -78.61367222189911,     -77.89151335098109,     -77.99854578005444,     -78.12499511405836,     -77.34585940736146,     -78.27762465340649,     -78.85710385889680,
			    -77.11546448332642,     -77.27420070195693,     -78.66417019326146,     -76.94535630969840,     -77.60525793048810,     -77.89644068440026,     -77.95009629561375,     -77.47451320893208,
			    -73.51512646579846,     -76.98031159471938,     -75.88140974114081,     -78.31160036127432,     -77.79164089234050,     -77.58599013089531,     -77.69000389648010,     -77.63384375481459,
			    -77.86542324847812,     -77.85828193114527,     -77.51595615652762,     -78.35614318223668,     -76.64013906444424,     -76.67941615864648,     -76.98477833172983,     -76.99379462976296,
			    -76.51898925046255,     -77.73034442364713,     -77.30188666626603,     -76.25031978709389,     -77.40252440820049,     -77.57689468044350,     -77.18615272577283,     -77.21123249887019,
			    -77.63937790834134,     -77.53070148215677,     -77.89020134603182,     -77.79120299478234,     -76.68778015789320,     -77.53924911481323,     -76.96669172754567,     -76.85781548613728,
			    -75.91489581912470,     -75.92281695611236,     -77.28945346184700,     -77.44019734495906,     -75.60854668633293,     -76.48216200467292,     -76.20897817840081,     -77.04210583981428,
			    -76.99317400330818,     -77.02549795173547,     -75.83900357663990,     -75.97062569099067,     -77.39780859836414,     -76.70313106195070,     -75.86926377499127,     -72.24506183135618,
			    -76.23851504277839,     -75.90345105314091,     -75.53456826746918,     -76.26858723648553,     -76.57363441142881,     -76.36439010064461,     -76.39687076903137,     -75.56430675981024,
			    -71.97985437881111,     -76.28326895732201,     -76.12795289017495,     -74.88148833978944,     -72.04422909393116,     -75.41439648759494,     -75.87472162518110,     -76.23319946533664,
			    -75.96315705620245,     -75.82337922012510,     -76.05849712945754,     -75.82220348743876,     -76.06468014192608,     -75.84885982447132,     -75.18686142918537,     -71.78019130293868,
			    -75.72111318168444,     -75.11062124193448,     -74.91111341690322,     -75.56813110614888,     -75.92840888638406,     -74.80621002623232,     -75.51811621454901,     -74.91918315140099,
			    -75.68063991169657,     -74.84686245038289,     -75.62192568740090,     -70.31024521321487,     -71.22206200624495,     -75.75950742230515,     -75.54140500708438,     -75.41698950569551,
			    -75.24602606513423,     -75.34181036840658,     -75.12125715588131,     -75.10487637836230,     -75.31768359603275,     -70.00501430628840,     -70.14355135986600,     -71.02157429000553,
			    -74.85177989672256,     -74.38139230623682,     -74.80037526612807,     -75.20005949407680,     -75.65007823182893,     -71.62809534644329,     -74.82600771311526,     -75.03162020299868,
			    -69.83331958751081,     -70.24688371971692,     -75.26703438650401,     -75.48516674042730,     -74.68750103469580,     -71.47749750158160,     -71.60707343115304,     -75.20716022870290,
			    -71.12508705697392,     -74.58876464909534,     -69.27169202142629,     -71.46807256559991,     -70.04164135150324,     -69.98350169811587,     -74.98753723356720,     -70.38491931355908,
			    -75.36874813111282,     -74.83654244019870,     -70.95896898633856,     -74.18469489204098,     -71.21380892677229,     -71.34234298702069,     -74.28070809176711,     -70.86198202002940,
			    -71.11742144437588,     -74.11903190686624,     -75.02875719705068,     -74.83466891018970,     -74.66929384092391,     -73.75127518935486,     -69.77877508593718,     -69.56335698554004,
			    -69.11254704320073,     -70.66955800090837,     -70.78433908287845,     -74.60209873734127,     -70.61122050484666,     -70.85391827131096,     -69.70508403081804,     -69.10571286753655,
			    -69.39988438439657,     -69.22321944776033,     -73.43177487157551,     -69.76264510156211,     -70.34823324507029,     -70.45476778554000,     -70.58410126244033,     -69.44178914987098,
			    -70.10233481383167,     -68.91278826409268,     -69.92336763894386,     -69.98001956767484,     -69.13723749620281,     -68.65494510532255,     -69.01880958500496,     -69.46640307354250,
			    -68.38245063847992,     -70.03229252587590,     -68.68080534724791,     -70.11708167800633,     -69.36720441284169,     -70.09417198886595,     -68.53128467942160,     -68.42262791218758,
			    -71.13773520890494,     -71.28357017901409,     -68.43250916165990,     -70.36234450498753,     -69.57848288888570,     -70.64730624365528,     -68.86256549425171,     -70.69410239227530,
			    -69.58801396093801,     -68.29631274866703,     -68.97947100782142,     -69.31038038696776,     -68.68258974430664,     -70.09834761675174,     -70.55746257558064,     -68.15855018186515,
			    -67.95665237825011,     -68.37724291541468,     -69.31407535699181,     -69.96186563973420,     -70.68339978913261,     -67.69889502798209,     -69.36439053240305,     -70.22778468201577,
			    -69.16318325488166,     -69.03141185767112,     -69.71396665333644,     -70.28347925992239,     -69.06703258669249,     -70.16312371859513,     -67.57527246464466,     -70.11012160252203,
			    -68.79829358720353,     -69.07725746580532,     -68.88691919808937,     -69.20232859228645,     -69.19317965865702,     -68.97000678976539,     -68.88311072161527,     -68.77822531641732,
			    -69.55147846861296,     -68.51690920606258,     -68.78008863416046,     -68.58519547642160,     -69.78267376612098,     -67.37297054605158,     -70.10887957236608,     -67.98934435388186,
			    -68.58878761025761,     -68.51749993299848,     -68.68261918678554,     -69.15014172112373,     -68.29887091620628,     -68.36896379318415,     -68.28723570141796,     -68.99476031016906,
			    -67.69429438262675,     -68.38458752687015,     -69.56593668861076,     -68.88327667250056,     -68.77390745853445,     -68.44845850396120,     -68.75062657802421,     -68.21900714983065,
			    -69.06359619332768,     -68.12889971814148,     -68.96864141568823,     -68.08734847686351,     -68.60217836856110,     -68.82173655702167,     -67.35105304539525,     -68.96460806581881,
			    -67.71469607474569,     -69.71534236509672,     -69.64817309841072,     -67.46327807283014,     -69.46104010956368,     -68.55179170981343,     -67.99052500196925,     -68.01573507774896,
			    -67.85148841419823,     -68.53295102586083,     -66.62593431495318,     -68.39639985433097,     -68.16819477078842,     -68.61360435526082,     -67.35526994400459,     -67.54984837283044,
			    -67.70023517610366,     -69.25768432401732,     -66.73390505870512,     -69.37938662536897,     -70.05426132218798,     -68.16738085905854,     -67.70779679403633,     -69.78267331494172,
			    -69.63563189081250,     -68.97658397241463,     -70.43030089034501,     -66.96823419456582,     -68.28566342914061,     -67.48235559016294,     -67.60554212758320,     -69.06010533797007,
			    -68.15359614844748,     -68.19906307255803,     -67.39245870403161,     -68.56740695831225,     -69.59583689638836,     -67.46167650679998,     -67.42638475197484,     -69.63912424484270,
			    -70.21574936561568,     -66.47497532472602,     -67.41244796827287,     -67.29732081289696,     -69.02707118825647,     -66.24309915545224,     -68.98702283054718,     -67.74436267204939,
			    -66.16899430489954,     -69.32460815646960,     -66.97380439005693,     -69.34797607129921,     -67.45580947632124,     -66.16813298286006,     -66.98370930528094,     -68.93198844190259,
			    -70.05966415946932,     -70.28776777019485,     -66.33177075178840,     -66.66830351001121,     -68.76835355105106,     -69.00445206242594,     -65.53823057346138,     -70.08302526891076,
			    -66.57381052288791,     -69.58827760212077,     -68.67698551667857,     -69.82907329938521,     -69.59451161612277,     -69.91170235964215,     -69.48693048304246,     -66.06287905103851,
			    -69.32992496096615,     -65.45256262904196,     -66.50934251598655,     -69.77623827238767,     -69.63902599973241,     -68.33629988627825,     -69.47972948929637,     -68.76804943362410,
			    -67.65953005683771,     -69.41011607949622,     -69.56619628060822,     -69.35270438146938,     -69.54762624580032,     -65.98431802112779,     -67.06828070827270,     -69.43351793377319,
			    -68.25310234026001,     -68.21863153989520,     -69.04353616980970,     -65.68206632116270,     -68.58640127340158,     -66.68529107899866,     -66.90340289930639,     -67.45275105923017,
			    -66.62568105626184,     -65.67987277103570,     -69.09229360914765,     -65.14413225836140,     -65.52728292770892,     -66.40197876276238,     -65.37777838814068,     -66.25494151421333,
			    -68.07084019398582,     -65.60461221437649,     -57.57973720829673,     -68.69382248316735,     -68.27129701617510,     -64.31547597116935,     -67.78535402236493,     -65.53769752699047,
			    -66.68672828682551,     -64.67447549409490,     -66.53917992866664,     -64.84078679592533,     -66.22609611612967,     -65.15685275506510,     -65.22255192836138,     -68.43088879820064,
			    -64.67635160746593,     -65.45099497967100,     -64.96116175347619,     -66.19915509338696,     -67.78195050187578,     -65.83985544990293,     -64.37130214128338,     -57.44325653908821,
			    -67.05091707098522,     -64.65899968750301,     -60.45474130135046,     -68.19990105460750,     -64.37273339328699,     -65.77255838585371,     -64.60434078029236,     -64.50440741614497,
			    -64.11037851180699,     -67.65799238507289,     -60.19057130343808,     -67.76384515487473,     -66.70070440263113,     -63.64081000106585,     -63.64288617583371,     -67.28269479523281,
			    -61.44463794695690,     -66.74776059281248,     -61.19085919531681,     -59.43196431607516,     -61.18244203672710,     -57.10438147035058,     -61.08905034644763,     -57.21878496663653,
			    -60.92882233340617,     -65.68199860113134,     -56.61528838811135,     -60.82656800451635,     -56.64548601953074,     -65.96498410027883,     -65.60694575778442,     -60.42128219834139,
			    -56.13774752550957,     -66.45934360001007,     -56.06771726944744,     -62.04304745199309,     -65.62530880984448,     -60.16759930499469,     -66.02470403837906,     -60.06599051752167,
			    -66.12017217376496,     -56.06998120053347,     -65.88295671508833,     -61.84329172489984,     -60.66998711721433,     -65.18243160619880,     -65.62953219002976,     -66.00098128953785,
			    -65.57881219023199,     -55.51909754607038,     -65.83286705959642,     -65.60396482248217,     -65.48850886552381,     -55.67741889993331,     -61.83686285226511,     -55.60589462764646,
			    -60.83970822226552,     -60.17406532325873,     -65.38999630571679,     -60.08133546004306,     -60.23629900211291,     -60.78436191010160,     -60.26525036713162,     -59.90782625445571,
			    -61.49715654292569,     -60.04354599314379,     -61.27615673662012,     -61.40483645420026,     -61.50047397441652,     -60.92404204967154,     -61.00455905106411,     -59.83998661243668,
			    -60.86694915466570,     -60.76098025301409,     -61.19705528040717,     -60.91146751300139,     -59.91212490872199,     -60.91529888029228,     -60.00426549518815,     -59.19533689604138,
			    -59.51513101716242,     -60.98895719669247,     -58.89244938557390,     -60.97374486062900,     -59.89705532417543,     -60.40949410663407,     -59.61438810152039,     -60.47430539591564,
			    -59.50376485582129,     -59.45800230243820,     -59.42983343937973,     -59.19881787880341,     -59.82064529891701,     -58.14041620005100,     -60.51855486348487,     -58.68536555502477,
			    -58.94841597749993,     -58.38119937650868,     -58.78635496415418,     -58.68880145213283,     -59.79939218493411,     -58.48074299815380,     -59.21693790279454,     -59.54257037723461,
			    -58.39221827468558,     -59.02328588095841,     -57.92505779589759,     -59.86657461721053,     -59.44447652835318,     -55.48489340901560,     -59.26133041502970,     -59.35965858251791,
			    -59.19453446542340,     -59.45821800712184,     -57.98292130824457,     -59.60381098449988,     -58.35079842297197,     -58.90798674421139,     -58.22024865124554,     -59.22921121602127,
			    -58.95679872834509,     -57.96421123039581,     -59.11328642303794,     -58.98236968890850,     -58.22334726847782,     -57.64226991582849,     -58.90918656720022,     -59.37283388923479,
			    -59.16734753296237,     -58.03349041680236,     -58.93629837892888,     -58.87029453539102,     -57.70178442535929,     -58.12167562687939,     -58.93688172728955,     -58.81782075402366,
			    -58.77334261144883,     -58.65151124012706,     -58.35005749274136,     -58.13380931838030,     -58.45294847470480,     -58.28242313608647,     -56.82520156531910,     -57.78098122689324,
			    -57.38541199832506,     -57.46575612262286,     -57.87693628120886,     -57.04498573320542,     -57.37113510054262,     -56.77886836814727,     -57.44247911179196,     -57.53792526949167,
			    -56.22929436297287,     -57.30097997325646,     -57.04684568452692,     -56.99667319216330,     -47.99263207295245,     -56.35468642906536,     -56.99492065519860,     -57.17097022122974,
			    -57.00271148917821,     -56.65859437792935,     -47.81064786487055,     -59.44405962040021,     -51.29071980527621,     -56.55997467825771,     -59.13919505804697,     -55.73970612166717
		};
		
		final double Epsilon = 1e-6;
		
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
				System.out.print(String.format("%22.14f", energy));
			}
		}
	}
	
	private static void benchmarkOps(SearchProblem search, SimpleConfSpace simpleConfSpace, List<ScoredConf> confs) {
		
		// settings
		final int[] numThreadsList = { 1, 2, 4, 8, 16 };//, 32, 64, 128 };
		final int[] numStreamsPerGpuList = { 1, 2, 4, 8, 16, 32, 64 };//, 128 };
		
		GpuStreamPool.printPoolSize = false;
		final int MaxNumGpus = Gpus.get().getGpus().size();
		
		Factory<ForcefieldInteractions,Molecule> interactionsFactory = (mol) -> FFInterGen.makeFullConf(search.confSpace, search.shellResidues, mol);
		ForcefieldParams ffparams = new ForcefieldParams();
		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.addTemplates(simpleConfSpace)
			.setParallelism(Parallelism.makeCpu(4))
			.build();
		ResPairCache resPairCache = new ResPairCache(ffparams, connectivity);
		
		for (int numThreads : numThreadsList) {
			ConfMinimizer minimizer = new CpuConfMinimizer.Builder(ffparams, interactionsFactory, search.confSpace)
				.setNumThreads(numThreads)
				.build();
			benchmarkOps(minimizer, confs, "CPU: " + numThreads);
			minimizer.cleanup();
		}
		
		// benchmark cpu
		for (int numThreads : numThreadsList) {
			new EnergyCalculator.Builder(simpleConfSpace, ffparams)
				.setType(EnergyCalculator.Type.Cpu)
				.setParallelism(Parallelism.makeCpu(numThreads))
				.setResPairCache(resPairCache)
				.use((fragEcalc) -> {
					benchmarkOps(simpleConfSpace, fragEcalc, confs, "CPU Residue: " + numThreads);
				});
		}
		
		// benchmark gpu
		for (int numGpus = 1; numGpus <= MaxNumGpus; numGpus*=2) {
			final int fNumGpus = numGpus;
		
			// benchmark cuda
			for (int numStreams : numStreamsPerGpuList) {
				ConfMinimizer minimizer = new GpuConfMinimizer.Builder(ffparams, interactionsFactory, search.confSpace)
					.setGpuInfo(GpuConfMinimizer.Type.Cuda, numGpus, numStreams)
					.build();
				benchmarkOps(minimizer, confs, "Cuda: " + numGpus + "x" + numStreams);
				minimizer.cleanup();
			}
			
			// benchmark cuda ccd
			for (int numStreams : numStreamsPerGpuList) {
				ConfMinimizer minimizer = new GpuConfMinimizer.Builder(ffparams, interactionsFactory, search.confSpace)
					.setGpuInfo(GpuConfMinimizer.Type.CudaCCD, numGpus, numStreams)
					.build();
				benchmarkOps(minimizer, confs, "Cuda CCD: " + numGpus + "x" + numStreams);
				minimizer.cleanup();
			}
			
			// benchmark residue cuda
			for (int numStreams : numStreamsPerGpuList) {
				new EnergyCalculator.Builder(simpleConfSpace, ffparams)
					.setType(EnergyCalculator.Type.ResidueCuda)
					.setParallelism(Parallelism.make(4, numGpus, numStreams))
					.setResPairCache(resPairCache)
					.use((fragEcalc) -> {
						benchmarkOps(simpleConfSpace, fragEcalc, confs, "Residue Cuda: " + fNumGpus + "x" + numStreams);
					});
			}
			
			// benchmark residue cuda ccd
			for (int numStreams : numStreamsPerGpuList) {
				new EnergyCalculator.Builder(simpleConfSpace, ffparams)
					.setType(EnergyCalculator.Type.ResidueCudaCCD)
					.setParallelism(Parallelism.make(4, numGpus, numStreams))
					.setResPairCache(resPairCache)
					.use((fragEcalc) -> {
						benchmarkOps(simpleConfSpace, fragEcalc, confs, "Residue Cuda CCD: " + fNumGpus + "x" + numStreams);
					});
			}
		}
	}

	private static void benchmarkOps(SimpleConfSpace confSpace, EnergyCalculator ecalc, List<ScoredConf> confs, String name) {
		
		System.out.print(String.format("\nBenchmarking %30s ...", name));
		
		final long MaxTimeNs = 10*TimeFormatter.NSpS;
		final int MaxNumConfs = confs.size();
		
		int p = ecalc.tasks.getParallelism();
		ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();
		
		// warmup
		for (int r=0; r<2; r++) {
			for (int i=0; i<p; i++) {
				confEcalc.calcEnergyAsync(confs.get(0), (energy) -> { /* don't care */ });
			}
		}
		
		// benchmark
		int numConfs = 0;
		Stopwatch stopwatch = new Stopwatch().start();
		while (numConfs < MaxNumConfs && stopwatch.getTimeNs() < MaxTimeNs) {
			for (int i=0; i<p; i++) {
				if (numConfs < MaxNumConfs) {
					confEcalc.calcEnergyAsync(confs.get(numConfs++), (energy) -> { /* don't care */ });
				}
			}
		}
		confEcalc.tasks.waitForFinish();
		stopwatch.stop();
		
		System.out.print(String.format(" confs: %5d, precise timing: %9s, ops: %7.2f",
			numConfs,
			stopwatch.getTime(TimeUnit.MILLISECONDS),
			numConfs/stopwatch.getTimeS()
		));
	}
	
	private static void benchmarkOps(ConfMinimizer minimizer, List<ScoredConf> confs, String name) {
		
		System.out.print(String.format("\nBenchmarking %30s ...", name));
		
		minimizer.setReportProgress(false);
		
		final long MaxTimeNs = 10*TimeFormatter.NSpS;
		final int MaxNumConfs = confs.size();
		
		ConfMinimizer.Async asyncMin = minimizer.getAsync();
		int p = asyncMin.getTasks().getParallelism();
		
		// warmup
		for (int r=0; r<2; r++) {
			for (int i=0; i<p; i++) {
				asyncMin.minimizeAsync(confs.get(0), (energy) -> { /* don't care */ });
			}
		}
		
		// benchmark
		int numConfs = 0;
		Stopwatch stopwatch = new Stopwatch().start();
		while (numConfs < MaxNumConfs && stopwatch.getTimeNs() < MaxTimeNs) {
			for (int i=0; i<p; i++) {
				if (numConfs < MaxNumConfs) {
					asyncMin.minimizeAsync(confs.get(numConfs++), (energy) -> { /* don't care */ });
				}
			}
		}
		asyncMin.getTasks().waitForFinish();
		stopwatch.stop();
		
		System.out.print(String.format(" confs: %5d, precise timing: %9s, ops: %7.2f",
			numConfs,
			stopwatch.getTime(TimeUnit.MILLISECONDS),
			numConfs/stopwatch.getTimeS()
		));
	}
}
