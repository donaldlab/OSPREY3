package edu.duke.cs.osprey.minimization;

import java.io.File;
import java.util.ArrayList;

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
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.cuda.ContextPool;
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
		for (int i=0; i<10; i++) {
			conf = tree.nextConf();
		}
		
		// get an arbitrary energy function
		ParameterizedMoleculeCopy pmol = new ParameterizedMoleculeCopy(search.confSpace);
		GpuEnergyFunctionGenerator egen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), new ContextPool(1));
		GpuForcefieldEnergy efunc = egen.fullConfEnergy(search.confSpace, search.shellResidues, pmol.getCopiedMolecule());
		MoleculeModifierAndScorer f = new MoleculeModifierAndScorer(efunc, search.confSpace, new RCTuple(conf.getAssignments()), pmol);
		
		System.out.println("dofs: " + f.getNumDOFs());
		
		// regular minimization
		System.out.println("minimizing CCD...");
		Minimizer.NeedsCleanup ccdMinimizer = new SimpleCCDMinimizer(f);
		//Minimizer ccdMinimizer = new CCDMinimizer(f, false);
		Stopwatch ccdStopwatch = new Stopwatch().start();
		ccdMinimizer.minimize();
		double ccdEnergy = efunc.getEnergy();
		System.out.println(String.format("finished in %s, energy: %12.6f", ccdStopwatch.stop().getTime(1), ccdEnergy));
		ccdMinimizer.cleanup();
		
		// experimental minimization
		System.out.println("minimizing quasi-Newton...");
		Minimizer.NeedsCleanup qnMinimizer = new QuasiNewtonMinimizer(f);
		Stopwatch qnStopwatch = new Stopwatch().start();
		qnMinimizer.minimize();
		double qnEnergy = efunc.getEnergy();
		System.out.println(String.format("finished in %s, energy: %12.6f", qnStopwatch.stop().getTime(1), qnEnergy));
		checkEnergy(ccdEnergy, qnEnergy);
		qnMinimizer.cleanup();
		
		f.cleanup();
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
