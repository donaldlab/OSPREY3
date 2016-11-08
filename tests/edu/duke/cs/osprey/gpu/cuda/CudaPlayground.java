package edu.duke.cs.osprey.gpu.cuda;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.ojalgo.optimisation.MathProgSysModel;

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
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.ForcefieldInteractionsGenerator;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.cuda.kernels.ForcefieldKernelCuda;
import edu.duke.cs.osprey.gpu.cuda.kernels.ForcefieldKernelOneBlockCuda;
import edu.duke.cs.osprey.gpu.cuda.kernels.SubForcefieldsKernelCuda;
import edu.duke.cs.osprey.minimization.LineSearcher;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;
import edu.duke.cs.osprey.minimization.SurfingLineSearcher;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.parallelism.TimingThread;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class CudaPlayground extends TestBase {
	
	public static void main(String[] args)
	throws Exception {
		
		initDefaultEnvironment();
		
		// NOTE: samples and such here:
		// https://github.com/jcuda/jcuda-samples/tree/master/JCudaSamples/src/main/java/jcuda
		
		// info on dynamic parallelism:
		// http://docs.nvidia.com/cuda/cuda-c-programming-guide/#cuda-dynamic-parallelism
		
		//forcefield();
		//linesearch();
		ccd();
		//subForcefields();
	}
	
	private static SearchProblem makeSearch()
	throws IOException {
		
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
		final boolean doBenchmarks = false;
		
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
		
		Stopwatch gpuOldStopwatch = null;
		if (doBenchmarks) {
			
			// benchmark
			// 1024 threads
			// 16 blocks: ~14.7k ops
			System.out.println("\nbenchmarking old GPU dof...");
			gpuOldStopwatch = new Stopwatch().start();
			for (int i=0; i<numRuns; i++) {
				cudaMof.getValForDOF(d, x.get(d));
			}
			System.out.println(String.format("finished in %s, %.1f ops, speedup: %.1fx\n",
				gpuOldStopwatch.stop().getTime(1),
				numRuns/gpuOldStopwatch.getTimeS(),
				(double)cpuStopwatch.getTimeNs()/gpuOldStopwatch.getTimeNs()
			));
		}
		
		// cleanup
		cudaMof.cleanup();
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
		ForcefieldInteractionsGenerator intergen = new ForcefieldInteractionsGenerator();
		ForcefieldInteractions interactions = intergen.makeFullConf(search.confSpace, search.shellResidues, cudaMol.getCopiedMolecule());
		BigForcefieldEnergy bigff = new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct);
		ForcefieldKernelOneBlockCuda ffkernel = new ForcefieldKernelOneBlockCuda(stream, bigff, dofs);
		
		// restore coords
		cudaMof.setDOFs(x);
		ffkernel.uploadCoordsAsync();
		
		// check accuracy
		checkEnergy(cpuEnergy, ffkernel.calcEnergySync());
		checkEnergy(cpuEnergyD, ffkernel.calcEnergyDofSync(d));
		
		double dihedralRadians = Math.toRadians(x.get(d));
		checkEnergy(cpuEnergyD, ffkernel.poseAndCalcEnergyDofSync(d, dihedralRadians));
		
		// restore coords
		cudaMof.setDOFs(x);
		ffkernel.uploadCoordsAsync();
		
		Stopwatch gpuOneBlockStopwatch = null;
		if (doBenchmarks) {
			
			// benchmark pose and dof energy
			// 1024 threads
			// 1 block: ~5.2k ops
			// 2 blocks: ~8.1k ops
			// 4 blocks: ~10.1k ops
			// 8 blocks: ~12.2k ops
			// 16 blocks: ~13.6k ops
			System.out.println("\nbenchmarking one-block GPU dof...");
			gpuOneBlockStopwatch = new Stopwatch().start();
			for (int i=0; i<numRuns; i++) {
				ffkernel.uploadCoordsAsync();
				ffkernel.calcEnergyDofSync(d);
			}
			System.out.println(String.format("finished in %s, %.1f ops, speedup over cpu: %.1fx, speedup over old: %.1fx\n",
				gpuOneBlockStopwatch.stop().getTime(1),
				numRuns/gpuOneBlockStopwatch.getTimeS(),
				(double)cpuStopwatch.getTimeNs()/gpuOneBlockStopwatch.getTimeNs(),
				(double)gpuOldStopwatch.getTimeNs()/gpuOneBlockStopwatch.getTimeNs()
			));
		}
		
		// cleanup
		cudaPool.release(stream);
		ffkernel.cleanup();
		cudaPool.cleanup();
		
		if (doBenchmarks) {
			
			class StreamThread extends TimingThread {
				
				private GpuEnergyFunctionGenerator egen;
				private ParameterizedMoleculeCopy mol;
				private GpuForcefieldEnergy efunc;
				private MoleculeModifierAndScorer mof;
				
				public StreamThread(int i, GpuEnergyFunctionGenerator egen) {
					super("stream-" + i);
					this.egen = egen;
				}

				@Override
				protected void warmup() {
					
					mol = new ParameterizedMoleculeCopy(search.confSpace);
					efunc = egen.fullConfEnergy(search.confSpace, search.shellResidues, mol.getCopiedMolecule());
					mof = new MoleculeModifierAndScorer(efunc, search.confSpace, tuple, mol);
					
					for (int i=0; i<100; i++) {
						mof.getValForDOF(d, x.get(d));
					}
				}

				@Override
				protected void time() {
					for (int i=0; i<numRuns; i++) {
						mof.getValForDOF(d, x.get(d));
					}
				}
				
				@Override
				protected void cleanup() {
					efunc.cleanup();
					mof.cleanup();
				}
			}
			
			int[] numStreamsList = { 1, 2 };//{ 1, 2, 4, 8, 16, 32, 64, 128 };
			
			long stream1TimeNs = 0;
			
			for (int numStreams : numStreamsList) {
				
				// benchmark streams
				// cpu: 2.2k
				
				// wait: spin - 1024 threads
				// 1 block:     1:5.2    2:9.9    4:16.9   8:16.2
				// 2 blocks:    1:8.5    2:15.8   4:23.5   8:20.1
				// 4 blocks:    1:10.8   2:20.1   4:27.0   8:23.0
				// 8 blocks:    1:13.2   2:24.1   4:27.8   8:23.9
				// 16 blocks:   1:13.4   2:26.3   4:27.8   8:24.1
				
				// wait: yield - 128 threads
				// 4 blocks:    1:6.6    2:11.3   4:20.4   8:23.9
				// 16 blocks:   1:11.6   2:19.3   4:28.0   8:25.8
				// n blocks:    1:17.1   2:28.0   4:33.2   8:26.2
				// wait: yield - 256 threads
				// 4 blocks:    1:9.8    2:16.3   4:27.6   8:25.7
				// 16 blocks:   1:13.1   2:22.3   4:27.9   8:25.3
				// n blocks:    1:16.8   2:30.1   4:28.5   8:25.2
				// wait: yield - 512 threads
				// 4 blocks:    1:10.1   2:18.8   4:27.3   8:25.4
				// 16 blocks:   1:14.3   2:25.3   4:28.2   8:25.6
				// n blocks:    1:16.9   2:31.1   4:29.1   8:26.1
				// wait: yield - 1024 threads
				// 1 block:     1:5.2    2:9.9    4:18.5   8:24.5   16:21.2
				// 2 blocks:    1:8.6    2:15.3   4:26.3   8:24.6
				// 4 blocks:    1:10.7   2:19.8   4:27.4   8:25.1
				// 8 blocks:    1:13.2   2:23.5   4:28.1   8:25.6
				// 16 blocks:   1:14.8   2:28.3   4:29.6   8:25.5
				// n blocks:    1:16.8   2:26.4   4:29.1   8:28.8
				
				// wait: yield - n blocks
				// 128 threads:    1:17.1   2:28.0   4:33.2   8:26.2
				// 256 threads:    1:16.8   2:30.1   4:28.5   8:25.2
				// 512 threads:    1:16.9   2:31.1   4:29.1   8:26.1
				// 1024 threads:   1:16.8   2:26.4   4:29.1   8:28.8
				
				// wait: blocking sync is really bad
				// yield < spin < blocking sync
				// more blocks => more better
				
				// sample to find tradeoff between threads and streams
				
				int numGpus = 1;
				GpuStreamPool streamPool = new GpuStreamPool(numGpus, divUp(numStreams, numGpus));
				GpuEnergyFunctionGenerator streamsEgen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), streamPool);
				List<TimingThread> threads = new ArrayList<>();
				for (int i=0; i<numStreams; i++) {
					threads.add(new StreamThread(i, streamsEgen));
				}
				System.out.println("benchmarking " + numStreams + " GPU streams...");
				Stopwatch gpuStreamsStopwatch = TimingThread.timeThreads(threads);
				if (numStreams == 1) {
					stream1TimeNs = gpuStreamsStopwatch.getTimeNs();
				}
				System.out.println(String.format("finished in %s, %.1f ops, speedup over cpu: %.1fx, speedup over 1 stream: %.1fx",
					gpuStreamsStopwatch.getTime(1),
					numStreams*numRuns/gpuStreamsStopwatch.getTimeS(),
					(double)numStreams*cpuStopwatch.getTimeNs()/gpuStreamsStopwatch.getTimeNs(),
					(double)numStreams*stream1TimeNs/gpuStreamsStopwatch.getTimeNs()
				));
				streamsEgen.cleanup();
			}
		}
	}
	
	private static void linesearch()
	throws IOException {
		
		SearchProblem search = makeSearch();
		RCTuple tuple = getConf(search, 0);
		ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
		
		final int NumRuns = 1000;
		final int d = 0;
		final boolean doBenchmarks = false;
		
		// init cpu side
		ParameterizedMoleculeCopy cpuMol = new ParameterizedMoleculeCopy(search.confSpace);
		EnergyFunction cpuEfunc = EnvironmentVars.curEFcnGenerator.fullConfEnergy(search.confSpace, search.shellResidues, cpuMol.getCopiedMolecule());
		MoleculeModifierAndScorer cpuMof = new MoleculeModifierAndScorer(cpuEfunc, search.confSpace, tuple, cpuMol);
		
		DoubleMatrix1D x = DoubleFactory1D.dense.make(cpuMof.getNumDOFs());
		ObjectiveFunction.DofBounds dofBounds = new ObjectiveFunction.DofBounds(cpuMof.getConstraints());
		dofBounds.getCenter(x);
		
		SurfingLineSearcher cpuLineSearcher = new SurfingLineSearcher();
		ObjectiveFunction.OneDof cpuFd = new ObjectiveFunction.OneDof(cpuMof, d);
		cpuLineSearcher.init(cpuFd);
		
		// restore coords
		cpuMof.setDOFs(x);
		
		double cpuXdstar = cpuLineSearcher.search(x.get(d));
		double cpuFxdstar = cpuFd.getValue(cpuXdstar);
		
		Stopwatch cpuStopwatch = null;
		if (doBenchmarks) {
			
			// benchmark cpu side: ~0.35k ops
			System.out.println("\nbenchmarking cpu...");
			cpuStopwatch = new Stopwatch().start();
			for (int i=0; i<NumRuns; i++) {
				cpuLineSearcher.search(x.get(d));
			}
			System.out.println(String.format("finished in %8s, ops: %5.0f\n",
				cpuStopwatch.stop().getTime(TimeUnit.MILLISECONDS),
				NumRuns/cpuStopwatch.getTimeS()
			));
		}
		
		/* TEMP
		// init opencl side
		ParameterizedMoleculeCopy openclMol = new ParameterizedMoleculeCopy(search.confSpace);
		GpuEnergyFunctionGenerator openclEgen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), new GpuQueuePool(1, 1));
		GpuForcefieldEnergy openclEfunc = openclEgen.fullConfEnergy(search.confSpace, search.shellResidues, openclMol.getCopiedMolecule());
		MoleculeModifierAndScorer openclMof = new MoleculeModifierAndScorer(openclEfunc, search.confSpace, tuple, openclMol);
		ObjectiveFunction.OneDof openclFd = new ObjectiveFunction.OneDof(openclMof, d);
		
		System.out.println("benchmarking opencl original...");
		SurfingLineSearcher openclOriginalLineSearcher = new SurfingLineSearcher();
		openclOriginalLineSearcher.init(openclFd);
		Stopwatch openclOriginalStopwatch = new Stopwatch().start();
		for (int i=0; i<NumRuns; i++) {
			xdstar = openclOriginalLineSearcher.search(xd);
		}
		System.out.println(String.format("finished in %8s, ops: %5.0f, speedup over cpu: %.1fx, dxd*: %.6f",
			openclOriginalStopwatch.stop().getTime(TimeUnit.MILLISECONDS),
			NumRuns/openclOriginalStopwatch.getTimeS(),
			(float)cpuStopwatch.getTimeNs()/openclOriginalStopwatch.getTimeNs(),
			xdstar - cpuXdstar
		));
		
		System.out.println("benchmarking opencl pipelined...");
		OpenCLSurfingLineSearcher openclPipelinedLineSearcher = new OpenCLSurfingLineSearcher();
		openclPipelinedLineSearcher.init(openclFd);
		Stopwatch openclPipelinedStopwatch = new Stopwatch().start();
		for (int i=0; i<NumRuns; i++) {
			xdstar = openclPipelinedLineSearcher.search(xd);
		}
		System.out.println(String.format("finished in %8s, ops: %5.0f, speedup over cpu: %.1fx, dxd*: %.6f",
			openclPipelinedStopwatch.stop().getTime(TimeUnit.MILLISECONDS),
			NumRuns/openclPipelinedStopwatch.getTimeS(),
			(float)cpuStopwatch.getTimeNs()/openclPipelinedStopwatch.getTimeNs(),
			xdstar - cpuXdstar
		));
		
		// cleanup
		openclMof.cleanup();
		openclEfunc.cleanup();
		openclEgen.cleanup();
		*/
		
		// init cuda side
		ParameterizedMoleculeCopy cudaMol = new ParameterizedMoleculeCopy(search.confSpace);
		GpuEnergyFunctionGenerator cudaEgen = new GpuEnergyFunctionGenerator(ffparams, new GpuStreamPool(1));
		GpuForcefieldEnergy cudaEfunc = cudaEgen.fullConfEnergy(search.confSpace, search.shellResidues, cudaMol.getCopiedMolecule());
		MoleculeModifierAndScorer cudaMof = new MoleculeModifierAndScorer(cudaEfunc, search.confSpace, tuple, cudaMol);
		ObjectiveFunction.OneDof cudaFd = new ObjectiveFunction.OneDof(cudaMof, d);
		
		SurfingLineSearcher cudaOriginalLineSearcher = new SurfingLineSearcher();
		cudaOriginalLineSearcher.init(cudaFd);
		
		// restore coords
		cudaMof.setDOFs(x);
		
		// check accuracy
		double xdstar = cudaOriginalLineSearcher.search(x.get(d));
		checkEnergy(cpuXdstar, xdstar);
		checkEnergy(cpuFxdstar, cudaFd.getValue(xdstar));
		
		Stopwatch cudaOriginalStopwatch = null;
		if (doBenchmarks) {
			
			// benchmark cuda side: ~2.4k ops
			System.out.println("\nbenchmarking cuda original...");
			cudaOriginalStopwatch = new Stopwatch().start();
			for (int i=0; i<NumRuns; i++) {
				cudaOriginalLineSearcher.search(x.get(d));
			}
			System.out.println(String.format("finished in %8s, ops: %5.0f, speedup over cpu: %.1fx\n",
				cudaOriginalStopwatch.stop().getTime(TimeUnit.MILLISECONDS),
				NumRuns/cudaOriginalStopwatch.getTimeS(),
				(float)cpuStopwatch.getTimeNs()/cudaOriginalStopwatch.getTimeNs()
			));
		}
		
		// cleanup
		cudaMof.cleanup();
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
		ForcefieldInteractionsGenerator intergen = new ForcefieldInteractionsGenerator();
		ForcefieldInteractions interactions = intergen.makeFullConf(search.confSpace, search.shellResidues, cudaMol.getCopiedMolecule());
		BigForcefieldEnergy bigff = new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct);
		ForcefieldKernelOneBlockCuda ffkernel = new ForcefieldKernelOneBlockCuda(stream, bigff, dofs);
		
		// restore coords
		cudaMof.setDOFs(x);
		ffkernel.uploadCoordsAsync();
		
		double xdRadians = Math.toRadians(x.get(d));
		double xdminRadians = Math.toRadians(dofBounds.getMin(d));
		double xdmaxRadians = Math.toRadians(dofBounds.getMax(d));
		double stepRadians = Math.toRadians(0.25);
		
		// check accuracy
		ForcefieldKernelOneBlockCuda.LinesearchResult linesearchResult = ffkernel.linesearchSync(d, xdRadians, xdminRadians, xdmaxRadians, stepRadians);
		checkEnergy(cpuXdstar, Math.toDegrees(linesearchResult.dihedralRadians));
		checkEnergy(cpuFxdstar, linesearchResult.energy);
		
		Stopwatch cudaOneBlockStopwatch = null;
		if (doBenchmarks) {
			
			// benchmark: ~1.1k ops
			System.out.println("\nbenchmarking cuda one-block...");
			cudaOneBlockStopwatch = new Stopwatch().start();
			for (int i=0; i<NumRuns; i++) {
				ffkernel.linesearchSync(d, xdRadians, xdminRadians, xdmaxRadians, stepRadians);	
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
		ffkernel.cleanup();
		cudaPool.cleanup();
		
		if (doBenchmarks) {
		
			class StreamThread extends TimingThread {
				
				private GpuStreamPool streams;
				private GpuStream stream;
				private ForcefieldKernelOneBlockCuda ffkernel;
				
				public StreamThread(int i, GpuStreamPool streams) {
					super("stream-" + i);
					this.streams = streams;
				}
				
				@Override
				protected void init()
				throws Exception {
					stream = streams.checkout();
					ForcefieldInteractionsGenerator intergen = new ForcefieldInteractionsGenerator();
					ForcefieldInteractions interactions = intergen.makeFullConf(search.confSpace, search.shellResidues, cudaMol.getCopiedMolecule());
					BigForcefieldEnergy bigff = new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct);
					ffkernel = new ForcefieldKernelOneBlockCuda(stream, bigff, dofs);
				}

				@Override
				protected void warmup() {
					for (int i=0; i<10; i++) {
						ffkernel.linesearchSync(d, xdRadians, xdminRadians, xdmaxRadians, stepRadians);	
					}
				}

				@Override
				protected void time() {
					for (int i=0; i<NumRuns; i++) {
						ffkernel.linesearchSync(d, xdRadians, xdminRadians, xdmaxRadians, stepRadians);	
					}
				}
				
				@Override
				protected void cleanup() {
					ffkernel.cleanup();
					streams.release(stream);
				}
			}
			
			int[] numStreamsList = { 2, 4, 8, 16, 32 };
			
			for (int numStreams : numStreamsList) {
				
				int numGpus = 1;
				cudaPool = new GpuStreamPool(numGpus, divUp(numStreams, numGpus));
				
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
		
		SimpleCCDMinimizer cpuMinimizer = new SimpleCCDMinimizer(cpuMof, new Factory<LineSearcher,Void>() {
			@Override
			public LineSearcher make(Void context) {
				return new SurfingLineSearcher();
			}
		});
		
		// restore coords
		cpuMof.setDOFs(x);
		
		DoubleMatrix1D cpuXstar = cpuMinimizer.minimize();
		double cpuEnergy = cpuMof.getValue(cpuXstar);
		
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
		
		SimpleCCDMinimizer cudaMinimizer = new SimpleCCDMinimizer(cudaMof, new Factory<LineSearcher,Void>() {
			@Override
			public LineSearcher make(Void context) {
				return new SurfingLineSearcher();
			}
		});
		
		// restore coords
		cudaMof.setDOFs(x);
		
		// check accuracy
		DoubleMatrix1D cudaXstar = cudaMinimizer.minimize();
		System.out.println(String.format("max xd dist: %8.6f", maxxddist(cpuXstar, cudaXstar, false)));
		checkEnergy(cpuEnergy, cpuMof.getValue(cpuXstar));
		
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
		cudaMof.cleanup();
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
		ForcefieldInteractionsGenerator intergen = new ForcefieldInteractionsGenerator();
		ForcefieldInteractions interactions = intergen.makeFullConf(search.confSpace, search.shellResidues, cudaMol.getCopiedMolecule());
		BigForcefieldEnergy bigff = new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct);
		ForcefieldKernelOneBlockCuda ffkernel = new ForcefieldKernelOneBlockCuda(stream, bigff, dofs);
		
		// restore coords
		cudaMof.setDOFs(x);
		ffkernel.uploadCoordsAsync();
		
		// check accuracy
		ForcefieldKernelOneBlockCuda.CCDResult result = ffkernel.ccdSync(x, dofBounds);
		System.out.println(String.format("max xd dist: %8.6f", maxxddist(cpuXstar, result.x, true)));
		checkEnergy(cpuEnergy, result.energy);
		
		Stopwatch cudaOneBlockStopwatch = null;
		if (doBenchmarks) {
			
			// benchmark: ~5 ops
			System.out.println("\nbenchmarking cuda one-block...");
			cudaOneBlockStopwatch = new Stopwatch().start();
			for (int i=0; i<NumRuns; i++) {
				ffkernel.uploadCoordsAsync();
				ffkernel.ccdSync(x, dofBounds);
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
		ffkernel.cleanup();
		cudaPool.cleanup();
		
		if (doBenchmarks) {
		
			class StreamThread extends TimingThread {
				
				private GpuStreamPool streams;
				private GpuStream stream;
				private ForcefieldKernelOneBlockCuda ffkernel;
				
				public StreamThread(int i, GpuStreamPool streams) {
					super("stream-" + i);
					this.streams = streams;
				}
				
				@Override
				protected void init()
				throws Exception {
					stream = streams.checkout();
					ForcefieldInteractionsGenerator intergen = new ForcefieldInteractionsGenerator();
					ForcefieldInteractions interactions = intergen.makeFullConf(search.confSpace, search.shellResidues, cudaMol.getCopiedMolecule());
					BigForcefieldEnergy bigff = new BigForcefieldEnergy(ffparams, interactions, BufferTools.Type.Direct);
					ffkernel = new ForcefieldKernelOneBlockCuda(stream, bigff, dofs);
				}

				@Override
				protected void warmup() {
					for (int i=0; i<2; i++) {
						ffkernel.uploadCoordsAsync();
						ffkernel.ccdSync(x, dofBounds);
					}
				}

				@Override
				protected void time() {
					for (int i=0; i<NumRuns; i++) {
						ffkernel.uploadCoordsAsync();
						ffkernel.ccdSync(x, dofBounds);
					}
				}
				
				@Override
				protected void cleanup() {
					ffkernel.cleanup();
					streams.release(stream);
				}
			}
			
			int[] numStreamsList = { 2, 4, 8, 16, 32 };
			
			for (int numStreams : numStreamsList) {
				
				int numGpus = 1;
				cudaPool = new GpuStreamPool(numGpus, divUp(numStreams, numGpus));
				
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
	
	private static double maxxddist(DoubleMatrix1D a, DoubleMatrix1D b, boolean convertToDegrees) {
		double maxDist = 0;
		assert (a.size() == b.size());
		int numDofs = a.size();
		for (int d=0; d<numDofs; d++) {
			double ad = a.get(d);
			double bd = b.get(d);
			if (convertToDegrees) {
				bd = Math.toDegrees(bd);
			}
			maxDist = Math.max(maxDist, Math.abs(ad - bd));
		}
		return maxDist;
	}

	private static void subForcefields()
	throws Exception {
		
		SearchProblem search = makeSearch();
		RCTuple tuple = getConf(search, 0);
		ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
		
		// init cpu side
		ParameterizedMoleculeCopy cpuMol = new ParameterizedMoleculeCopy(search.confSpace);
		EnergyFunction cpuEfunc = EnvironmentVars.curEFcnGenerator.fullConfEnergy(search.confSpace, search.shellResidues, cpuMol.getCopiedMolecule());
		MoleculeModifierAndScorer cpuMof = new MoleculeModifierAndScorer(cpuEfunc, search.confSpace, tuple, cpuMol);
		
		DoubleMatrix1D x = DoubleFactory1D.dense.make(cpuMof.getNumDOFs());
		ObjectiveFunction.DofBounds dofBounds = new ObjectiveFunction.DofBounds(cpuMof.getConstraints());
		
		dofBounds.getCenter(x);
		double cpuEnergy = cpuMof.getValue(x);
		
		// init cuda side
		ParameterizedMoleculeCopy cudaMol = new ParameterizedMoleculeCopy(search.confSpace);
		GpuEnergyFunctionGenerator cudaEgen = new GpuEnergyFunctionGenerator(ffparams, new GpuStreamPool(1));
		GpuForcefieldEnergy cudaEfunc = cudaEgen.fullConfEnergy(search.confSpace, search.shellResidues, cudaMol.getCopiedMolecule());
		MoleculeModifierAndScorer cudaMof = new MoleculeModifierAndScorer(cudaEfunc, search.confSpace, tuple, cudaMol);
		
		// init the forcefield kernel
		checkEnergy(cpuEnergy, cudaMof.getValue(x));
		
		List<FreeDihedral> dofs = new ArrayList<>();
		for (DegreeOfFreedom dof : cudaMof.getDOFs()) {
			dofs.add((FreeDihedral)dof);
		}
		ForcefieldKernelCuda cudaFFKernel = (ForcefieldKernelCuda)cudaEfunc.getKernel();
		SubForcefieldsKernelCuda subForcefieldsKernel = new SubForcefieldsKernelCuda(cudaFFKernel, dofs);
		
		subForcefieldsKernel.calcEnergies(x, 0.25);
		checkEnergy(cpuEnergy, cpuEnergy + subForcefieldsKernel.getFxdmOffset(0));
		checkEnergy(cpuEnergy, cpuEnergy + subForcefieldsKernel.getFxdpOffset(0));
		
		// cleanup
		cudaMof.cleanup();
		cudaEfunc.cleanup();
		cudaEgen.cleanup();
	}

	private static void checkEnergy(double cpuEnergy, double gpuEnergy) {
		System.out.println(String.format("cpu: %12.6f   gpu: %12.6f   err: %.12f", cpuEnergy, gpuEnergy, Math.abs(cpuEnergy - gpuEnergy)));
	}
	
	private static int divUp(int a, int b) {
		return (a + b - 1)/b;
	}
}
