package edu.duke.cs.osprey.minimization;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.ojalgo.matrix.decomposition.LU;
import org.ojalgo.matrix.store.PrimitiveDenseStore;

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
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.cuda.ContextPool;
import edu.duke.cs.osprey.gpu.cuda.kernels.ForcefieldKernelCuda;
import edu.duke.cs.osprey.gpu.cuda.kernels.SubForcefieldsKernelCuda;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class MinimizationPlayground extends TestBase {
	
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
		
		// get the ith conformation
		search.pruneMat = new PruningMatrix(search.confSpace, 1000);
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 4, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		ScoredConf conf = null;
		for (int i=0; i<100000; i++) {
			conf = tree.nextConf();
		}
		// NOTE: conf 100 apparently has better energies in the QN minimizer??? =D
		
		// get an arbitrary energy function
		ParameterizedMoleculeCopy pmol = new ParameterizedMoleculeCopy(search.confSpace);
		GpuEnergyFunctionGenerator egen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), new ContextPool(1));
		GpuForcefieldEnergy efunc = egen.fullConfEnergy(search.confSpace, search.shellResidues, pmol.getCopiedMolecule());
		MoleculeModifierAndScorer f = new MoleculeModifierAndScorer(efunc, search.confSpace, new RCTuple(conf.getAssignments()), pmol);
		
		System.out.println("dofs: " + f.getNumDOFs());
		
		// warm up the efunc
		{
			DoubleMatrix1D x = DoubleFactory1D.dense.make(f.getNumDOFs());
			new ObjectiveFunction.DofBounds(f.getConstraints()).getCenter(x);
			for (int i=0; i<10; i++) {
				f.getValue(x);
			}
		}
		
		// original ccd minimization
		System.out.println("minimizing CCD...");
		CCDMinimizer ccdMinimizer = new CCDMinimizer(f, false);
		Stopwatch ccdStopwatch = new Stopwatch().start();
		ccdMinimizer.minimize();
		double ccdEnergy = efunc.getEnergy();
		System.out.println(String.format("finished in %s, energy: %12.6f", ccdStopwatch.stop().getTime(1), ccdEnergy));
		
		// simple ccd minimization
		System.out.println("minimizing simple CCD...");
		SimpleCCDMinimizer simpleCcdMinimizer = new SimpleCCDMinimizer(f);
		Stopwatch simpleCcdStopwatch = new Stopwatch().start();
		simpleCcdMinimizer.minimize();
		double simpleCcdEnergy = efunc.getEnergy();
		System.out.println(String.format("finished in %s, energy: %12.6f", simpleCcdStopwatch.stop().getTime(1), simpleCcdEnergy));
		checkEnergy(ccdEnergy, simpleCcdEnergy);
		simpleCcdMinimizer.cleanup();
		
		// warm up the LU solver so it doesn't affect our timings
		LU<Double> solver = LU.PRIMITIVE.make();
		int m = f.getNumDOFs()*2 + 1;
		PrimitiveDenseStore A = PrimitiveDenseStore.FACTORY.makeEye(m, m);
		PrimitiveDenseStore b = PrimitiveDenseStore.FACTORY.makeZero(m, 1);
		solver.decompose(A);
		assert (solver.isSolvable());
		solver.solve(b);
		
		// quasi-newton minimization
		System.out.println("minimizing quasi-Newton...");
		QuasiNewtonMinimizer qnMinimizer = new QuasiNewtonMinimizer(f);
		Stopwatch qnStopwatch = new Stopwatch().start();
		qnMinimizer.minimize();
		double qnEnergy = efunc.getEnergy();
		System.out.println(String.format("finished in %s, energy: %12.6f", qnStopwatch.stop().getTime(1), qnEnergy));
		checkEnergy(ccdEnergy, qnEnergy);
		f.cleanup();
		
		// init cuda minimization directly
		ParameterizedMoleculeCopy cudaMol = new ParameterizedMoleculeCopy(search.confSpace);
		GpuEnergyFunctionGenerator cudaEgen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), new ContextPool(1));
		GpuForcefieldEnergy cudaEfunc = cudaEgen.fullConfEnergy(search.confSpace, search.shellResidues, cudaMol.getCopiedMolecule());
		MoleculeModifierAndScorer cudaMof = new MoleculeModifierAndScorer(cudaEfunc, search.confSpace, new RCTuple(conf.getAssignments()), cudaMol);
		List<FreeDihedral> dofs = new ArrayList<>();
		for (DegreeOfFreedom dof : cudaMof.getDOFs()) {
			dofs.add((FreeDihedral)dof);
		}
		ForcefieldKernelCuda cudaFFKernel = (ForcefieldKernelCuda)cudaEfunc.getKernel();
		SubForcefieldsKernelCuda subForcefieldsKernel = new SubForcefieldsKernelCuda(cudaFFKernel, dofs);
		QuasiNewtonMinimizer gpuQnMinimizer = new QuasiNewtonMinimizer(cudaMof, subForcefieldsKernel);
		
		// gpu quasi-newton minimization
		System.out.println("minimizing GPU quasi-Newton...");
		Stopwatch gpuQnStopwatch = new Stopwatch().start();
		gpuQnMinimizer.minimize();
		double gpuQnEnergy = cudaEfunc.getEnergy();
		System.out.println(String.format("finished in %s, energy: %12.6f", gpuQnStopwatch.stop().getTime(1), gpuQnEnergy));
		checkEnergy(ccdEnergy, gpuQnEnergy);
		
		// cleanup
		subForcefieldsKernel.cleanup();
		cudaEfunc.cleanup();
	}
	
	private static void checkEnergy(double expected, double observed) {
		
		final double Epsilon = 1e-6;
			
		double absErr = observed - expected;
		if (absErr > Epsilon) {
			
			System.out.println(String.format("low precision energy: exp:%12.8f  obs:%12.8f       absErr:%12.8f",
				expected, observed, absErr
			));
			
		} else if (absErr < -Epsilon) {
		
			System.out.println(String.format("improved energy: exp:%12.8f  obs:%12.8f  improvement:%12.8f",
				expected, observed, -absErr
			));
		}
	}
}
