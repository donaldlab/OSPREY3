package edu.duke.cs.osprey.gpu.cuda;

import java.io.File;
import java.io.IOException;
import java.nio.DoubleBuffer;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.cuda.kernels.TestDPKernel;
import edu.duke.cs.osprey.gpu.cuda.kernels.TestKernel;
import edu.duke.cs.osprey.minimization.CudaSurfingLineSearcher;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.minimization.SurfingLineSearcher;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
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
		
		// init CUDA
		Context context = new Context(Gpus.get().getGpus().get(0));
		
		dynamicParallelism(context);
		linesearch(context);
	
		context.cleanup();
	}
	
	private static void dynamicParallelism(Context context)
	throws IOException {
		
		final int NumElements = 122000;
		final int NumRuns = 1000;
	
		System.out.println("dynamic parallelism test: " + NumElements + " threads");
		
		// init host-loop kernel
		TestKernel hostKernel = new TestKernel(context, NumElements);
		DoubleBuffer a = hostKernel.getA();
		DoubleBuffer b = hostKernel.getB();
		a.clear();
		b.clear();
		for (int i=0; i<NumElements; i++) {
			a.put((double)i);
			b.put((double)i);
		}
		
		System.out.println("host-side loop...");
		Stopwatch hostStopwatch = new Stopwatch().start();
		for (int i=0; i<NumRuns; i++) {
			hostKernel.uploadAsync();
			hostKernel.runAsync();
			hostKernel.downloadSync();
		}
		System.out.println(String.format("finished in %8s, ops: %5.0f",
			hostStopwatch.stop().getTime(TimeUnit.MILLISECONDS),
			NumRuns/hostStopwatch.getTimeS()
		));
		hostKernel.cleanup();
		
		// init device-loop kernel
		TestDPKernel deviceKernel = new TestDPKernel(context, NumElements);
		a = deviceKernel.getA();
		b = deviceKernel.getB();
		a.clear();
		b.clear();
		for (int i=0; i<NumElements; i++) {
			a.put((double)i);
			b.put((double)i);
		}
		
		System.out.println("device-side loop...");
		Stopwatch deviceStopwatch = new Stopwatch().start();
		deviceKernel.uploadAsync();
		deviceKernel.runAsync();
		deviceKernel.downloadSync();
		System.out.println(String.format("finished in %8s, ops: %5.0f, speedup: %.1fx",
			deviceStopwatch.stop().getTime(TimeUnit.MILLISECONDS),
			NumRuns/deviceStopwatch.getTimeS(),
			(float)hostStopwatch.getTimeNs()/deviceStopwatch.getTimeNs()
		));
		deviceKernel.cleanup();
	}
	
	private static void linesearch(Context context)
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
		
		RCTuple tuple = new RCTuple();
		for (int i=0; i<search.confSpace.numPos; i++) {
			tuple = tuple.addRC(i, 0);
		}
		
		int d = 0;
		
		// init cpu side
		ParameterizedMoleculeCopy cpuMol = new ParameterizedMoleculeCopy(search.confSpace);
		EnergyFunction cpuEfunc = EnvironmentVars.curEFcnGenerator.fullConfEnergy(search.confSpace, search.shellResidues, cpuMol.getCopiedMolecule());
		MoleculeModifierAndScorer cpuMof = new MoleculeModifierAndScorer(cpuEfunc, search.confSpace, tuple, cpuMol);
		ObjectiveFunction.OneDof cpuFd = new ObjectiveFunction.OneDof(cpuMof, d);
		
		final int NumRuns = 1000;
		double xd = (cpuFd.getXMin() + cpuFd.getXMax())/2;
		double xdstar = 0;
		
		// benchmark cpu side
		System.out.println("benchmarking cpu...");
		SurfingLineSearcher cpuLineSearcher = new SurfingLineSearcher();
		cpuLineSearcher.init(cpuFd);
		Stopwatch cpuStopwatch = new Stopwatch().start();
		for (int i=0; i<NumRuns; i++) {
			xdstar = cpuLineSearcher.search(xd);
		}
		System.out.println(String.format("finished in %8s, ops: %5.0f",
			cpuStopwatch.stop().getTime(TimeUnit.MILLISECONDS),
			NumRuns/cpuStopwatch.getTimeS()
		));
		
		double cpuXdstar = xdstar;
		
		
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
		GpuEnergyFunctionGenerator cudaEgen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), new ContextPool(1));
		GpuForcefieldEnergy cudaEfunc = cudaEgen.fullConfEnergy(search.confSpace, search.shellResidues, cudaMol.getCopiedMolecule());
		MoleculeModifierAndScorer cudaMof = new MoleculeModifierAndScorer(cudaEfunc, search.confSpace, tuple, cudaMol);
		ObjectiveFunction.OneDof cudaFd = new ObjectiveFunction.OneDof(cudaMof, d);
		
		System.out.println("Num atom pairs: " + cudaEfunc.getSubset().getNumAtomPairs());
		
		// benchmark cuda side
		System.out.println("benchmarking cuda original...");
		SurfingLineSearcher cudaOriginalLineSearcher = new SurfingLineSearcher();
		cudaOriginalLineSearcher.init(cudaFd);
		Stopwatch cudaOriginalStopwatch = new Stopwatch().start();
		for (int i=0; i<NumRuns; i++) {
			xdstar = cudaOriginalLineSearcher.search(xd);
		}
		System.out.println(String.format("finished in %8s, ops: %5.0f, speedup over cpu: %.1fx, dxd*: %.6f",
			cudaOriginalStopwatch.stop().getTime(TimeUnit.MILLISECONDS),
			NumRuns/cudaOriginalStopwatch.getTimeS(),
			(float)cpuStopwatch.getTimeNs()/cudaOriginalStopwatch.getTimeNs(),
			xdstar - cpuXdstar
		));
		
		// benchmark cuda side
		System.out.println("benchmarking cuda pipelined...");
		CudaSurfingLineSearcher cudaPipelinedLineSearcher = new CudaSurfingLineSearcher();
		cudaPipelinedLineSearcher.init(cudaFd);
		Stopwatch cudaPipelinedStopwatch = new Stopwatch().start();
		for (int i=0; i<NumRuns; i++) {
			xdstar = cudaPipelinedLineSearcher.search(xd);
		}
		System.out.println(String.format("finished in %8s, ops: %5.0f, speedup over cpu: %.1fx, dxd*: %.6f",
			cudaPipelinedStopwatch.stop().getTime(TimeUnit.MILLISECONDS),
			NumRuns/cudaPipelinedStopwatch.getTimeS(),
			(float)cpuStopwatch.getTimeNs()/cudaPipelinedStopwatch.getTimeNs(),
			xdstar - cpuXdstar
		));
		
		// cleanup
		cudaMof.cleanup();
		cudaEfunc.cleanup();
		cudaEgen.cleanup();
	}
}
