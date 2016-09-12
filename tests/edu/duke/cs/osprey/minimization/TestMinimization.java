package edu.duke.cs.osprey.minimization;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

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
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.GpuQueuePool;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class TestMinimization extends TestBase {
	
	private static double[] ExpectedEnergies = {
		-107.01471465433335, -107.14427781940432, -106.79145713231975, -106.92139365967053, -106.28769308885211, -107.11801397703762,
		-106.67892206113300, -106.41908247351522, -106.89600279606412, -106.74468003314176, -106.45689550906734, -106.95533592350961
	};
	
	// NOTE: minimization energies aren't super precise
	private static final double Epsilon = 1e-6;
	
	private static SearchProblem search;
	private static Factory<EnergyFunction,Molecule> efuncgen;
	private static Factory<GpuForcefieldEnergy,Molecule> gpuefuncgen;
	private static List<ScoredConf> confs;
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
		
		// for these small problems, more than one thread is actually slower
		MultiTermEnergyFunction.setNumThreads(1);
		
		String aaNames = "ALA";
		String mutRes = "39 43";
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
		
		search = new SearchProblem(
			"test", "test/1CC8/1CC8.ss.pdb", 
			flexResList, allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null, false
		);
		
		// make the cpu energy stuff
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		efuncgen = new Factory<EnergyFunction,Molecule>() {
			@Override
			public EnergyFunction make(Molecule mol) {
				return egen.fullConfEnergy(search.confSpace, search.shellResidues, mol);
			}
		};
		
		// make the gpu energy stuff
		GpuEnergyFunctionGenerator gpuegen = new GpuEnergyFunctionGenerator(makeDefaultFFParams(), new GpuQueuePool(1, 2));
		gpuefuncgen = new Factory<GpuForcefieldEnergy,Molecule>() {
			@Override
			public GpuForcefieldEnergy make(Molecule mol) {
				return gpuegen.fullConfEnergy(search.confSpace, search.shellResidues, mol);
			}
		};
		
		// compute the energy matrix and pruning matrix
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, search.confSpace, search.shellResidues);
		search.emat = new SimpleEnergyMatrixCalculator(ecalc).calcEnergyMatrix(); 
		search.pruneMat = new PruningMatrix(search.confSpace, 1000);
		
		// build A* tree
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		// get the confs
		confs = new ArrayList<>();
		for (int i=0; i<ExpectedEnergies.length; i++) {
			confs.add(tree.nextConf());
		}
	}
	
	private void assertEnergies(List<EnergiedConf> econfs) {
		
		assertThat(econfs.size(), is(confs.size()));
		
		for (int i=0; i<confs.size(); i++) {
			
			ScoredConf conf = confs.get(i);
			EnergiedConf econf = econfs.get(i);
			
			assertThat(econf.getAssignments(), is(conf.getAssignments()));
			assertThat(econf.getScore(), is(conf.getScore()));
			assertThat(econf.getEnergy(), isRelatively(ExpectedEnergies[i], Epsilon));
		}
	}
	
	@Test
	public void testOldWay() {
		
		for (int i=0; i<ExpectedEnergies.length; i++) {
			double energy = search.minimizedEnergy(confs.get(i).getAssignments());
			assertThat(energy, isRelatively(ExpectedEnergies[i], Epsilon));
		}
	}
	
	@Test
	public void testMainThread() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		EnergyFunction efunc = efuncgen.make(search.confSpace.m);
		
		// minimize on main thread
		List<EnergiedConf> econfs = new ArrayList<>();
		for (ScoredConf conf : confs) {
			econfs.add(minimizer.minimize(search.confSpace.m, conf, efunc, search.confSpace));
		}
		
		assertEnergies(econfs);
	}
	
	@Test
	public void testMainThreadBatch() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		EnergyFunction efunc = efuncgen.make(search.confSpace.m);
		
		// minimize on main thread, in batch mode
		List<EnergiedConf> econfs = minimizer.minimize(search.confSpace.m, confs, efunc, search.confSpace);
		
		assertEnergies(econfs);
	}
	
	@Test
	public void testMainThreadTasks() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		
		// minimize on main thread, in batch mode
		List<EnergiedConf> econfs = minimizer.minimize(confs, efuncgen, search.confSpace);
		
		assertEnergies(econfs);
	}
	
	@Test
	public void testTaskThread() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(1);
		
		// minimize on main thread, in batch mode
		List<EnergiedConf> econfs = minimizer.minimize(confs, efuncgen, search.confSpace, tasks);
		
		tasks.stop();
		
		assertEnergies(econfs);
	}
		
	@Test
	public void test2TaskThreads() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(2);
		
		// minimize on main thread, in batch mode
		List<EnergiedConf> econfs = minimizer.minimize(confs, efuncgen, search.confSpace, tasks);
		
		tasks.stop();
		
		assertEnergies(econfs);
	}
	
	@Test
	public void testMainGpu() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		GpuForcefieldEnergy efunc = gpuefuncgen.make(search.confSpace.m);
		
		// minimize on main thread
		List<EnergiedConf> econfs = new ArrayList<>();
		for (ScoredConf conf : confs) {
			econfs.add(minimizer.minimize(search.confSpace.m, conf, efunc, search.confSpace));
		}
		
		efunc.cleanup();
		
		assertEnergies(econfs);
	}
	
	@Test
	public void testTaskGpu() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(1);
		
		// minimize on main thread, in batch mode
		List<EnergiedConf> econfs = minimizer.minimize(confs, gpuefuncgen, search.confSpace, tasks);
		
		tasks.stop();
		
		assertEnergies(econfs);
	}
	
	@Test
	public void test2TaskGpus() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(2);
		
		// minimize on main thread, in batch mode
		List<EnergiedConf> econfs = minimizer.minimize(confs, gpuefuncgen, search.confSpace, tasks);
		
		tasks.stop();
		
		assertEnergies(econfs);
	}
}