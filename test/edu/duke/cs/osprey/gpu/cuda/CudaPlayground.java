package edu.duke.cs.osprey.gpu.cuda;

import java.io.File;
import java.io.IOException;
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
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.FFInterGen;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MinimizingFragmentEnergyCalculator;
import edu.duke.cs.osprey.energy.MinimizingFragmentEnergyCalculator.Type;
import edu.duke.cs.osprey.energy.ResInterGen;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ResPairCache;
import edu.duke.cs.osprey.energy.forcefield.ResidueForcefieldEnergy;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.cuda.kernels.CCDKernelCuda;
import edu.duke.cs.osprey.gpu.cuda.kernels.ResidueForcefieldEnergyCuda;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TimingThread;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tupexp.LUTESettings;

@SuppressWarnings("unused")
public class CudaPlayground extends TestBase {
	
	// TODO: split into accuracy unit tests and benchmarks
	
	public static void main(String[] args)
	throws Exception {
		
		initDefaultEnvironment();
		
		// NOTE: samples and such here:
		// https://github.com/jcuda/jcuda-samples/tree/master/JCudaSamples/src/main/java/jcuda
		
		// info on dynamic parallelism:
		// http://docs.nvidia.com/cuda/cuda-c-programming-guide/#cuda-dynamic-parallelism
		
		//forcefield();
		//ccd();
		residueCcd();
	}
	
	private static SearchProblem makeSearch()
	throws IOException {
		
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
			search.emat = new SimpleEnergyMatrixCalculator.Cpu(2, EnvironmentVars.curEFcnGenerator.ffParams, search.confSpace, search.shellResidues).calcEnergyMatrix();
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
		}
		
		return search;
	}
	
	private static RCTuple getConf(SearchProblem search, int i) {
		
		// get the ith conformation
		search.pruneMat = new PruningMatrix(search.confSpace, 1000);
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 4, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		ScoredConf conf = null;
		i++;
		for (int j=0; j<i; j++) {
			conf = tree.nextConf();
		}
		return new RCTuple(conf.getAssignments());
	}
	
	private static void forcefield()
	throws Exception {
		
		SearchProblem search = makeSearch();
		RCTuple tuple = getConf(search, 0);
		ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
		
		final int numRuns = 10000;
		final int d = 0;
		final boolean doBenchmarks = true;
		
		// init cpu side
		ParameterizedMoleculeCopy cpuMol = new ParameterizedMoleculeCopy(search.confSpace);
		EnergyFunction cpuEfunc = EnvironmentVars.curEFcnGenerator.fullConfEnergy(search.confSpace, search.shellResidues, cpuMol.getCopiedMolecule());
		MoleculeModifierAndScorer cpuMof = new MoleculeModifierAndScorer(cpuEfunc, search.confSpace, tuple, cpuMol);
		
		DoubleMatrix1D x = DoubleFactory1D.dense.make(cpuMof.getNumDOFs());
		ObjectiveFunction.DofBounds dofBounds = new ObjectiveFunction.DofBounds(cpuMof.getConstraints());
		
		dofBounds.getCenter(x);
		double cpuEnergy = cpuMof.getValue(x);
		double cpuEnergyD = cpuMof.getValForDOF(d, x.get(d));
		
		Stopwatch cpuStopwatch = null;
		if (doBenchmarks) {
			
			// benchmark
			System.out.println("\nbenchmarking CPU...");
			cpuStopwatch = new Stopwatch().start();
			for (int i=0; i<numRuns; i++) {
				cpuMof.getValForDOF(d, x.get(d));
			}
			System.out.println(String.format("finished in %s, %.1f ops\n",
				cpuStopwatch.stop().getTime(1),
				numRuns/cpuStopwatch.getTimeS()
			));
		}
		
		System.out.println();
		
		// init cuda side
		ParameterizedMoleculeCopy cudaMol = new ParameterizedMoleculeCopy(search.confSpace);
		GpuEnergyFunctionGenerator cudaEgen = new GpuEnergyFunctionGenerator(ffparams, new GpuStreamPool(1));
		GpuForcefieldEnergy cudaEfunc = cudaEgen.fullConfEnergy(search.confSpace, search.shellResidues, cudaMol.getCopiedMolecule());
		MoleculeModifierAndScorer cudaMof = new MoleculeModifierAndScorer(cudaEfunc, search.confSpace, tuple, cudaMol);
		
		// full efunc
		System.out.println("atom pairs: " + cudaEfunc.getKernel().getSubset().getNumAtomPairs());
		double gpuEnergy = cudaMof.getValue(x);
		checkEnergy(cpuEnergy, gpuEnergy);
		
		// one dof
		System.out.println("d " + d + " atom pairs: " + ((GpuForcefieldEnergy)cudaMof.getEfunc(d)).getSubset().getNumAtomPairs());
		double gpuEnergyD = cudaMof.getValForDOF(d, x.get(d));
		checkEnergy(cpuEnergyD, gpuEnergyD);
		
		Stopwatch gpuStopwatch = null;
		if (doBenchmarks) {
			
			// benchmark
			// 1024 threads
			// 16 blocks: ~14.7k ops
			System.out.println("\nbenchmarking GPU dof...");
			gpuStopwatch = new Stopwatch().start();
			for (int i=0; i<numRuns; i++) {
				cudaMof.getValForDOF(d, x.get(d));
			}
			System.out.println(String.format("finished in %s, %.1f ops, speedup: %.1fx\n",
				gpuStopwatch.stop().getTime(1),
				numRuns/gpuStopwatch.getTimeS(),
				(double)cpuStopwatch.getTimeNs()/gpuStopwatch.getTimeNs()
			));
		}
		
		// cleanup
		cudaEfunc.cleanup();
		cudaEgen.cleanup();
	}

	private static void ccd()
	throws IOException {
		
		SearchProblem search = makeSearch();
		RCTuple tuple = getConf(search, 0);
		ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
		
		final int NumRuns = 10;
		final boolean doBenchmarks = true;
		
		// init cpu side
		ParameterizedMoleculeCopy cpuMol = new ParameterizedMoleculeCopy(search.confSpace);
		EnergyFunction cpuEfunc = EnvironmentVars.curEFcnGenerator.fullConfEnergy(search.confSpace, search.shellResidues, cpuMol.getCopiedMolecule());
		MoleculeModifierAndScorer cpuMof = new MoleculeModifierAndScorer(cpuEfunc, search.confSpace, tuple, cpuMol);
		
		DoubleMatrix1D x = DoubleFactory1D.dense.make(cpuMof.getNumDOFs());
		ObjectiveFunction.DofBounds dofBounds = new ObjectiveFunction.DofBounds(cpuMof.getConstraints());
		dofBounds.getCenter(x);
		
		SimpleCCDMinimizer cpuMinimizer = new SimpleCCDMinimizer();
		cpuMinimizer.init(cpuMof);
		
		// restore coords
		cpuMof.setDOFs(x);
		
		Minimizer.Result cpuResult = cpuMinimizer.minimize();
		
		Stopwatch cpuStopwatch = null;
		if (doBenchmarks) {
			
			// benchmark cpu side: ~2 ops
			System.out.println("\nbenchmarking cpu...");
			cpuStopwatch = new Stopwatch().start();
			for (int i=0; i<NumRuns; i++) {
				cpuMof.setDOFs(x);
				cpuMinimizer.minimize();
			}
			System.out.println(String.format("finished in %8s, ops: %5.0f\n",
				cpuStopwatch.stop().getTime(TimeUnit.MILLISECONDS),
				NumRuns/cpuStopwatch.getTimeS()
			));
		}
		
		// init cuda side
		ParameterizedMoleculeCopy cudaMol = new ParameterizedMoleculeCopy(search.confSpace);
		GpuEnergyFunctionGenerator cudaEgen = new GpuEnergyFunctionGenerator(ffparams, new GpuStreamPool(1));
		GpuForcefieldEnergy cudaEfunc = cudaEgen.fullConfEnergy(search.confSpace, search.shellResidues, cudaMol.getCopiedMolecule());
		MoleculeModifierAndScorer cudaMof = new MoleculeModifierAndScorer(cudaEfunc, search.confSpace, tuple, cudaMol);
		
		SimpleCCDMinimizer cudaMinimizer = new SimpleCCDMinimizer();
		cudaMinimizer.init(cudaMof);
		
		// restore coords
		cudaMof.setDOFs(x);
		
		// check accuracy
		Minimizer.Result cudaResult = cudaMinimizer.minimize();
		System.out.println(String.format("max xd dist: %8.6f", maxxddist(cpuResult.dofValues, cudaResult.dofValues)));
		checkEnergy(cpuResult.energy, cpuMof.getValue(cpuResult.dofValues));
		
		Stopwatch cudaOriginalStopwatch = null;
		if (doBenchmarks) {
			
			// benchmark cuda side: ~15 ops
			System.out.println("\nbenchmarking cuda original...");
			cudaOriginalStopwatch = new Stopwatch().start();
			for (int i=0; i<NumRuns; i++) {
				cudaMof.setDOFs(x);
				cudaMinimizer.minimize();
			}
			System.out.println(String.format("finished in %8s, ops: %5.0f, speedup over cpu: %.1fx\n",
				cudaOriginalStopwatch.stop().getTime(TimeUnit.MILLISECONDS),
				NumRuns/cudaOriginalStopwatch.getTimeS(),
				(float)cpuStopwatch.getTimeNs()/cudaOriginalStopwatch.getTimeNs()
			));
		}
		
		// cleanup
		cudaEfunc.cleanup();
		cudaEgen.cleanup();
		
		// collects the dofs
		List<FreeDihedral> dofs = new ArrayList<>();
		for (DegreeOfFreedom dof : cudaMof.getDOFs()) {
			dofs.add((FreeDihedral)dof);
		}
		
		// make the kernel directly
		GpuStreamPool cudaPool = new GpuStreamPool(1);
		GpuStream stream = cudaPool.checkout();
		ForcefieldInteractions interactions = FFInterGen.makeFullConf(search.confSpace, search.shellResidues, cudaMol.getCopiedMolecule());
		BigForcefieldEnergy bigff = new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct);
		MoleculeModifierAndScorer bigMof = new MoleculeModifierAndScorer(bigff, search.confSpace, tuple, cudaMol);
		CCDKernelCuda kernel = new CCDKernelCuda(stream);
		kernel.init(bigMof);
		
		// restore coords
		cudaMof.setDOFs(x);
		kernel.uploadCoordsAsync();
		
		// check accuracy
		Minimizer.Result result = kernel.runSync(x, dofBounds);
		System.out.println(String.format("max xd dist: %8.6f", maxxddist(cpuResult.dofValues, result.dofValues)));
		checkEnergy(cpuResult.energy, result.energy);
		
		Stopwatch cudaOneBlockStopwatch = null;
		if (doBenchmarks) {
			
			// benchmark: ~5 ops
			System.out.println("\nbenchmarking cuda one-block...");
			cudaOneBlockStopwatch = new Stopwatch().start();
			for (int i=0; i<NumRuns; i++) {
				kernel.uploadCoordsAsync();
				kernel.runSync(x, dofBounds);
			}
			System.out.println(String.format("finished in %8s, ops: %5.0f, speedup over cpu: %.1fx, speedup over original: %.1fx\n",
				cudaOneBlockStopwatch.stop().getTime(TimeUnit.MILLISECONDS),
				NumRuns/cudaOneBlockStopwatch.getTimeS(),
				(float)cpuStopwatch.getTimeNs()/cudaOneBlockStopwatch.getTimeNs(),
				(float)cudaOriginalStopwatch.getTimeNs()/cudaOneBlockStopwatch.getTimeNs()
			));
		}
		
		// cleanup
		cudaPool.release(stream);
		kernel.cleanup();
		cudaPool.cleanup();
		
		if (doBenchmarks) {
		
			class StreamThread extends TimingThread {
				
				private GpuStreamPool streams;
				private GpuStream stream;
				private CCDKernelCuda ffkernel;
				
				public StreamThread(int i, GpuStreamPool streams) {
					super("stream-" + i);
					this.streams = streams;
				}
				
				@Override
				protected void init()
				throws Exception {
					stream = streams.checkout();
					ForcefieldInteractions interactions = FFInterGen.makeFullConf(search.confSpace, search.shellResidues, cudaMol.getCopiedMolecule());
					BigForcefieldEnergy bigff = new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct);
					MoleculeModifierAndScorer bigMof = new MoleculeModifierAndScorer(bigff, search.confSpace, tuple, cudaMol);
					ffkernel = new CCDKernelCuda(stream);
					ffkernel.init(bigMof);
				}

				@Override
				protected void warmup() {
					for (int i=0; i<2; i++) {
						ffkernel.uploadCoordsAsync();
						ffkernel.runSync(x, dofBounds);
					}
				}

				@Override
				protected void time() {
					for (int i=0; i<NumRuns; i++) {
						ffkernel.uploadCoordsAsync();
						ffkernel.runSync(x, dofBounds);
					}
				}
				
				@Override
				protected void cleanup() {
					ffkernel.cleanup();
					streams.release(stream);
				}
			}
			
			int[] numStreamsList = { 1, 4, 8, 12, 16, 20, 24, 28, 32 };
			// 512 threads
			// 16 streams: 75.5, 82.3, 70.1 ops
			// 20 streams: 67.5, 81.2, 60.4 ops
			// 24 streams: 83.0, 61.8, 74.6 ops
			// 28 streams: 63.7, 74.8, 78.0 ops
			// 32 streams: 74.7, 82.2, 78.3 ops
			// wow, lots of volatility... not clear what's best
			
			for (int numStreams : numStreamsList) {
				
				int numGpus = 1;
				cudaPool = new GpuStreamPool(numGpus, MathTools.divUp(numStreams, numGpus));
				
				List<TimingThread> threads = new ArrayList<>();
				for (int i=0; i<numStreams; i++) {
					threads.add(new StreamThread(i, cudaPool));
				}
				System.out.println("benchmarking " + numStreams + " GPU streams...");
				Stopwatch gpuStreamsStopwatch = TimingThread.timeThreads(threads);
				System.out.println(String.format("finished in %s, %.1f ops, speedup over cpu: %.1fx, speedup over 1 stream: %.1fx",
					gpuStreamsStopwatch.getTime(1),
					numStreams*NumRuns/gpuStreamsStopwatch.getTimeS(),
					(double)numStreams*cpuStopwatch.getTimeNs()/gpuStreamsStopwatch.getTimeNs(),
					(double)numStreams*cudaOneBlockStopwatch.getTimeNs()/gpuStreamsStopwatch.getTimeNs()
				));
				
				cudaPool.cleanup();
			}
		}
	}
	
	private static void residueCcd()
	throws IOException {
		
		SearchProblem search = makeSearch();
		RCTuple tuple = getConf(search, 0);
		ForcefieldParams ffparams = new ForcefieldParams();
		
		// also make a simple conf space
		Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
		strand.flexibility.get(39).setLibraryRotamers("ALA").setContinuous();
		strand.flexibility.get(43).setLibraryRotamers("ALA").setContinuous();
		strand.flexibility.get(40).setLibraryRotamers().setContinuous();
		strand.flexibility.get(41).setLibraryRotamers().setContinuous();
		strand.flexibility.get(42).setLibraryRotamers().setContinuous();
		strand.flexibility.get(44).setLibraryRotamers().setContinuous();
		strand.flexibility.get(45).setLibraryRotamers().setContinuous();
		SimpleConfSpace simpleConfSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		assertConfSpacesMatch(search.confSpace, simpleConfSpace);
		
		// minimize on CPU side
		ParameterizedMoleculeCopy cpuMol = new ParameterizedMoleculeCopy(search.confSpace);
		EnergyFunction cpuEfunc = EnvironmentVars.curEFcnGenerator.fullConfEnergy(search.confSpace, search.shellResidues, cpuMol.getCopiedMolecule());
		MoleculeModifierAndScorer cpuMof = new MoleculeModifierAndScorer(cpuEfunc, search.confSpace, tuple, cpuMol);
		DoubleMatrix1D x = DoubleFactory1D.dense.make(cpuMof.getNumDOFs());
		ObjectiveFunction.DofBounds dofBounds = new ObjectiveFunction.DofBounds(cpuMof.getConstraints());
		dofBounds.getCenter(x);
		SimpleCCDMinimizer cpuMinimizer = new SimpleCCDMinimizer();
		cpuMinimizer.init(cpuMof);
		
		// restore coords
		cpuMof.setDOFs(x);
		
		Minimizer.Result cpuResult = cpuMinimizer.minimize();
		System.out.println(String.format("CPU energy: %12.6f", cpuResult.energy));
		
		// minimize with residue Cuda CCD
		MinimizingFragmentEnergyCalculator gpuFragEcalc = new MinimizingFragmentEnergyCalculator.Builder(simpleConfSpace, ffparams)
			.setType(Type.ResidueCudaCCD)
			.setParallelism(Parallelism.makeGpu(1, 1))
			.build();
		ResidueInteractions inters = new EnergyPartition.Traditional().makeFragment(simpleConfSpace, null, tuple);
		double gpuEnergy = gpuFragEcalc.calcEnergy(tuple, inters);
		
		// TEMP: try some simpler minimizations
		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.setConfSpace(simpleConfSpace)
			.setParallelism(Parallelism.makeCpu(4))
			.build();
		ResPairCache resPairCache = new ResPairCache(ffparams, connectivity);
		GpuStreamPool streams = new GpuStreamPool(1, 1);
		
		System.out.println();
		tuple = new RCTuple(2, 0, 3, 0);
		inters = ResInterGen.of(simpleConfSpace)
			.addIntra(2)
			.make();
		runCpuEfunc(simpleConfSpace, resPairCache, tuple, inters);
		runGpuEfunc(simpleConfSpace, resPairCache, tuple, inters, streams);
		gpuFragEcalc.calcEnergy(tuple, inters);
		
		System.out.println();
		tuple = new RCTuple(2, 0, 3, 0);
		inters = ResInterGen.of(simpleConfSpace)
			.addIntra(3)
			.make();
		runCpuEfunc(simpleConfSpace, resPairCache, tuple, inters);
		runGpuEfunc(simpleConfSpace, resPairCache, tuple, inters, streams);
		gpuFragEcalc.calcEnergy(tuple, inters);
		
		System.out.println();
		tuple = new RCTuple(2, 0, 3, 0);
		inters = ResInterGen.of(simpleConfSpace)
			.addIntra(2)
			.addInter(2, 3)
			.make();
		runCpuEfunc(simpleConfSpace, resPairCache, tuple, inters);
		runGpuEfunc(simpleConfSpace, resPairCache, tuple, inters, streams);
		gpuFragEcalc.calcEnergy(tuple, inters);
		
		System.out.println();
		tuple = new RCTuple(2, 0, 3, 0);
		inters = ResInterGen.of(simpleConfSpace)
			.addIntra(3)
			.addInter(2, 3)
			.make();
		runCpuEfunc(simpleConfSpace, resPairCache, tuple, inters);
		runGpuEfunc(simpleConfSpace, resPairCache, tuple, inters, streams);
		gpuFragEcalc.calcEnergy(tuple, inters);
		
		System.out.println();
		tuple = new RCTuple(2, 0, 3, 0);
		inters = ResInterGen.of(simpleConfSpace)
			.addInter(2, 3)
			.make();
		runCpuEfunc(simpleConfSpace, resPairCache, tuple, inters);
		runGpuEfunc(simpleConfSpace, resPairCache, tuple, inters, streams);
		gpuFragEcalc.calcEnergy(tuple, inters);
		
		System.out.println();
		tuple = new RCTuple(2, 0, 3, 0);
		inters = ResInterGen.of(simpleConfSpace)
			.addIntra(2)
			.addIntra(3)
			.addInter(2, 3)
			.make();
		runCpuEfunc(simpleConfSpace, resPairCache, tuple, inters);
		runGpuEfunc(simpleConfSpace, resPairCache, tuple, inters, streams);
		gpuFragEcalc.calcEnergy(tuple, inters);
		
		streams.cleanup();
	}
	
	private static void runCpuEfunc(SimpleConfSpace confSpace, ResPairCache resPairCache, RCTuple frag, ResidueInteractions inters) {
		ParametricMolecule pmol = confSpace.makeMolecule(frag);
		ObjectiveFunction.DofBounds bounds = confSpace.makeBounds(frag);
		ResidueForcefieldEnergy efunc = new ResidueForcefieldEnergy(resPairCache, inters, pmol.mol);
		MoleculeObjectiveFunction mof = new MoleculeObjectiveFunction(pmol, bounds, efunc);
		DoubleMatrix1D x = DoubleFactory1D.dense.make(bounds.size());
		bounds.getCenter(x);
		mof.setDOFs(x);
		System.out.println(String.format("CPU e:   %12.6f", efunc.getEnergy()));
	}
	
	private static void runGpuEfunc(SimpleConfSpace confSpace, ResPairCache resPairCache, RCTuple frag, ResidueInteractions inters, GpuStreamPool streams)
	throws IOException {
		ParametricMolecule pmol = confSpace.makeMolecule(frag);
		ObjectiveFunction.DofBounds bounds = confSpace.makeBounds(frag);
		ResidueForcefieldEnergyCuda efunc = new ResidueForcefieldEnergyCuda(streams, resPairCache, inters, pmol.mol);
		MoleculeObjectiveFunction mof = new MoleculeObjectiveFunction(pmol, bounds, efunc);
		DoubleMatrix1D x = DoubleFactory1D.dense.make(bounds.size());
		bounds.getCenter(x);
		mof.setDOFs(x);
		System.out.println(String.format("GPU e:   %12.6f", efunc.getEnergy()));
		efunc.cleanup();
	}
	
	private static double maxxddist(DoubleMatrix1D a, DoubleMatrix1D b) {
		double maxDist = 0;
		assert (a.size() == b.size());
		int numDofs = a.size();
		for (int d=0; d<numDofs; d++) {
			maxDist = Math.max(maxDist, Math.abs(a.get(d) - b.get(d)));
		}
		return maxDist;
	}

	private static void checkEnergy(double cpuEnergy, double gpuEnergy) {
		System.out.println(String.format("cpu: %12.6f   gpu: %12.6f   err: %.12f", cpuEnergy, gpuEnergy, Math.abs(cpuEnergy - gpuEnergy)));
	}
}
