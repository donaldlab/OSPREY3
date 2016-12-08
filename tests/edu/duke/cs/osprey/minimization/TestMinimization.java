package edu.duke.cs.osprey.minimization;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.ArrayList;
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
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.ForcefieldInteractionsGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class TestMinimization extends TestBase {
	
	private static double[] ExpectedEnergies = {
		-89.40966969379109,     -89.10792031500127,     -89.80959784194695,     -88.63999143548550,
		-89.12813398454155,     -89.50404412354314,     -88.39619842051209,     -88.88944810225344,
		-88.91539256575626,     -88.37401748235720,     -88.72521745741045,     -88.95852827540257,
		-88.56492542985106,     -89.13542390896973,     -88.39342805731060,     -88.61512935924652
	};
	
	private static final double Epsilon = 1e-8;
	
	private static SearchProblem search;
	private static ForcefieldParams ffparams;
	private static Factory<EnergyFunction,Molecule> efuncgen;
	private static Factory<GpuForcefieldEnergy,Molecule> openclEfuncgen;
	private static Factory<GpuForcefieldEnergy,Molecule> cudaEfuncgen;
	private static Factory<ForcefieldInteractions,Molecule> intergen;
	private static List<ScoredConf> confs;
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
		
		// for these small problems, more than one thread is actually slower
		MultiTermEnergyFunction.setNumThreads(1);
		
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
		
		search = new SearchProblem(
			"test", "test/1CC8/1CC8.ss.pdb", 
			resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null,
			false, new ArrayList<>()
		);
		
		ffparams = makeDefaultFFParams();
		
		// make energy function factories
		ForcefieldInteractionsGenerator ffintergen = new ForcefieldInteractionsGenerator();
		intergen = (mol) -> ffintergen.makeFullConf(search.confSpace, search.shellResidues, mol);
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		efuncgen = (mol) -> egen.interactionEnergy(intergen.make(mol));
		openclEfuncgen = (mol) -> new GpuForcefieldEnergy(ffparams, intergen.make(mol), new GpuQueuePool(1, 2));
		cudaEfuncgen = (mol) -> new GpuForcefieldEnergy(ffparams, intergen.make(mol), new GpuStreamPool(1, 2));
		
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
		final int numConfs = 16;
		confs = new ArrayList<>();
		for (int i=0; i<numConfs; i++) {
			confs.add(tree.nextConf());
		}
	}
	
	public static void main(String[] args) {
		
		before();
		
		// compute the expected energies
		for (int i=0; i<confs.size(); i++) {
			
			// get the objective function
			ParameterizedMoleculeCopy pmol = new ParameterizedMoleculeCopy(search.confSpace);
			EnergyFunction efunc = efuncgen.make(pmol.getCopiedMolecule());
			RCTuple tuple = new RCTuple(confs.get(i).getAssignments());
			MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, search.confSpace, tuple, pmol);
			
			// use the original CCD minimizer
			Minimizer.Result result = new CCDMinimizer(mof, false).minimize();
			
			// print the expected energy
			if (i > 0) {
				System.out.print(",");
			}
			System.out.print(i % 4 == 0 ? "\n" : " ");
			System.out.print(String.format("%22.14f", result.energy));
		}
	}
	
	private void assertEnergies(List<EnergiedConf> econfs) {
		
		assertThat(econfs.size(), is(confs.size()));
		
		for (int i=0; i<confs.size(); i++) {
			
			ScoredConf conf = confs.get(i);
			EnergiedConf econf = econfs.get(i);
			
			assertThat(econf.getAssignments(), is(conf.getAssignments()));
			assertThat(econf.getScore(), is(conf.getScore()));
			
			// penalize large errors, but not lower energies
			double absErr = econf.getEnergy() - ExpectedEnergies[i];
			assertThat(absErr, lessThanOrEqualTo(Epsilon));
		}
	}
	
	@Test
	public void testSearchProblemMinimizer() {
		
		for (int i=0; i<ExpectedEnergies.length; i++) {
			double energy = search.minimizedEnergy(confs.get(i).getAssignments());
			assertThat(energy, isAbsolutely(ExpectedEnergies[i], Epsilon));
		}
	}
	
	@Test
	public void testMainThread() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		EnergyFunction efunc = efuncgen.make(search.confSpace.m);
		
		// minimize on main thread
		List<EnergiedConf> econfs = new ArrayList<>();
		for (ScoredConf conf : confs) {
			econfs.add(minimizer.minimize(ParameterizedMoleculeCopy.makeNoCopy(search.confSpace), conf, efunc, search.confSpace));
		}
		
		assertEnergies(econfs);
	}
	
	@Test
	public void testMainThreadBatch() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		EnergyFunction efunc = efuncgen.make(search.confSpace.m);
		
		// minimize on main thread, in batch mode
		List<EnergiedConf> econfs = minimizer.minimize(ParameterizedMoleculeCopy.makeNoCopy(search.confSpace), confs, efunc, search.confSpace);
		
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
	public void testMainOpenCL() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		GpuForcefieldEnergy efunc = openclEfuncgen.make(search.confSpace.m);
		
		// minimize on main thread
		List<EnergiedConf> econfs = new ArrayList<>();
		for (ScoredConf conf : confs) {
			econfs.add(minimizer.minimize(ParameterizedMoleculeCopy.makeNoCopy(search.confSpace), conf, efunc, search.confSpace));
		}
		
		efunc.cleanup();
		
		assertEnergies(econfs);
	}
	
	@Test
	public void testTaskOpenCL() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(1);
		
		// minimize on main thread, in batch mode
		List<EnergiedConf> econfs = minimizer.minimize(confs, openclEfuncgen, search.confSpace, tasks);
		
		tasks.stop();
		
		assertEnergies(econfs);
	}
	
	@Test
	public void test2TaskOpenCL() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(2);
		
		// minimize on main thread, in batch mode
		List<EnergiedConf> econfs = minimizer.minimize(confs, openclEfuncgen, search.confSpace, tasks);
		
		tasks.stop();
		
		assertEnergies(econfs);
	}
	
	@Test
	public void testMainCuda() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		GpuForcefieldEnergy efunc = cudaEfuncgen.make(search.confSpace.m);
		
		// minimize on main thread
		List<EnergiedConf> econfs = new ArrayList<>();
		for (ScoredConf conf : confs) {
			econfs.add(minimizer.minimize(ParameterizedMoleculeCopy.makeNoCopy(search.confSpace), conf, efunc, search.confSpace));
		}
		
		efunc.cleanup();
		
		assertEnergies(econfs);
	}
	
	@Test
	public void testTaskCuda() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(1);
		
		// minimize on main thread, in batch mode
		List<EnergiedConf> econfs = minimizer.minimize(confs, cudaEfuncgen, search.confSpace, tasks);
		
		tasks.stop();
		
		assertEnergies(econfs);
	}
	
	@Test
	public void test2TaskCuda() {
		
		ConfMinimizer minimizer = new ConfMinimizer();
		
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(2);
		
		// minimize on main thread, in batch mode
		List<EnergiedConf> econfs = minimizer.minimize(confs, cudaEfuncgen, search.confSpace, tasks);
		
		tasks.stop();
		
		assertEnergies(econfs);
	}

	@Test
	public void testCpuConfMinimizer1Thread() {
		check(new CpuConfMinimizer(1, ffparams, intergen, search.confSpace));
	}
	
	@Test
	public void testCpuConfMinimizer2Threads() {
		check(new CpuConfMinimizer(2, ffparams, intergen, search.confSpace));
	}
	
	@Test
	public void testCudaConfMinmizer1Stream() {
		check(new GpuConfMinimizer(GpuConfMinimizer.Type.Cuda, 1, 1, ffparams, intergen, search.confSpace));
	}
	
	@Test
	public void testCudaConfMinmizer2Streams() {
		check(new GpuConfMinimizer(GpuConfMinimizer.Type.Cuda, 1, 2, ffparams, intergen, search.confSpace));
	}
	
	@Test
	public void testOpenCLConfMinmizer1Stream() {
		check(new GpuConfMinimizer(GpuConfMinimizer.Type.OpenCL, 1, 1, ffparams, intergen, search.confSpace));
	}
	
	@Test
	public void testOpenCLConfMinmizer2Streams() {
		check(new GpuConfMinimizer(GpuConfMinimizer.Type.OpenCL, 1, 2, ffparams, intergen, search.confSpace));
	}
	
	private void check(SpecializedConfMinimizer minimizer) {
		List<EnergiedConf> econfs = minimizer.minimize(confs);
		minimizer.cleanup();
		assertEnergies(econfs);
	}
}