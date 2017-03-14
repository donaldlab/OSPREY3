package edu.duke.cs.osprey.minimization;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.ForcefieldInteractionsGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class TestMinimization extends TestBase {
	
	private static Map<Boolean,Info> Infos;
	
	private static final double Epsilon = 1e-7;
	
	private static class Info {
		
		public ForcefieldParams ffparams;
		public SearchProblem search;
		public EnergyFunctionGenerator efuncgen;
		public Factory<ForcefieldInteractions,Molecule> intergen;
		public List<ScoredConf> confs;
		public double[] expectedEnergies;
		
		public Info(boolean doSolv, double[] expectedEnergies) {
			
			ffparams = makeDefaultFFParams();
			ffparams.doSolvationE = doSolv;
			
			this.expectedEnergies = expectedEnergies;
			
			efuncgen = new EnergyFunctionGenerator(ffparams, Double.POSITIVE_INFINITY, false);
			EnvironmentVars.curEFcnGenerator = efuncgen;
			
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
			
			// make energy function factories
			ForcefieldInteractionsGenerator ffintergen = new ForcefieldInteractionsGenerator();
			intergen = (mol) -> ffintergen.makeFullConf(search.confSpace, search.shellResidues, mol);
			
			// compute the energy matrix and pruning matrix
			search.emat = new SimpleEnergyMatrixCalculator.Cpu(2, ffparams, search.confSpace, search.shellResidues).calcEnergyMatrix();
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
	}
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
		
		// for these small problems, more than one thread is actually slower
		MultiTermEnergyFunction.setNumThreads(1);
		
		Infos = new HashMap<>();
		Infos.put(true, new Info(true, new double[] {
			-89.40966969380000,     -89.10792031501593,     -89.80959784194780,     -88.63999143525771,
			-89.12813398454430,     -89.50404412354364,     -88.39619842044286,     -88.88944810225273,
			-88.91539256575295,     -88.37401748235216,     -88.72521745744406,     -88.95852827535202,
			-88.56492542985737,     -89.13542390896987,     -88.39342805735203,     -88.61512935924370
		}));
		Infos.put(false, new Info(false, new double[] {
			-70.47295175891054,     -70.18816916335444,     -70.29217193373283,     -70.48982613389194,
			-69.42023183578743,     -69.91475960627933,     -70.35460308647797,     -69.19541718369999,
			-69.11087623590001,     -69.43654754385510,     -69.48055945486242,     -70.12045267097507,
			-69.42165499907402,     -68.97467845375986,     -69.34083168223796,     -69.74355354866111
		}));
	}
	
	public static void main(String[] args) {
		
		before();
		
		for (boolean doSolv : Arrays.asList(true, false)) {
		
			Info info = Infos.get(doSolv);
			
			System.out.println("Solvation? " + doSolv);
			
			// compute the expected energies
			for (int i=0; i<info.confs.size(); i++) {
				
				// get the objective function
				ParameterizedMoleculeCopy pmol = new ParameterizedMoleculeCopy(info.search.confSpace);
				EnergyFunction efunc = info.efuncgen.interactionEnergy(info.intergen.make(pmol.getCopiedMolecule()));
				RCTuple tuple = new RCTuple(info.confs.get(i).getAssignments());
				MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, info.search.confSpace, tuple, pmol);
				
				// use the original CCD minimizer
				Minimizer.Result result = new CCDMinimizer(mof, false).minimize();
				
				// print the expected energy
				if (i > 0) {
					System.out.print(",");
				}
				System.out.print(i % 4 == 0 ? "\n" : " ");
				System.out.print(String.format("%22.14f", result.energy));
			}
			
			System.out.println();
		}
	}
	
	@Test
	public void testSearchProblemMinimizer() {
		
		for (boolean doSolv : Arrays.asList(true, false)) {
			Info info = Infos.get(doSolv);
			
			for (int i=0; i<info.expectedEnergies.length; i++) {
				double energy = info.search.minimizedEnergy(info.confs.get(i).getAssignments());
				assertThat(energy, isAbsolutely(info.expectedEnergies[i], Epsilon));
			}
		}
	}
	
	@Test
	public void testCpuConfMinimizer1Thread() {
		check((ffparams, intergen, confSpace) -> new CpuConfMinimizer.Builder(ffparams, intergen, confSpace).build());
	}
	
	@Test
	public void testCpuConfMinimizer2Threads() {
		check((ffparams, intergen, confSpace) -> new CpuConfMinimizer.Builder(ffparams, intergen, confSpace).setNumThreads(2).build());
	}
	
	@Test
	public void testCudaConfMinmizer1Stream() {
		check((ffparams, intergen, confSpace) -> new GpuConfMinimizer.Builder(ffparams, intergen, confSpace).setGpuInfo(GpuConfMinimizer.Type.Cuda, 1, 1).build());
	}
	
	@Test
	public void testCudaConfMinmizer2Streams() {
		check((ffparams, intergen, confSpace) -> new GpuConfMinimizer.Builder(ffparams, intergen, confSpace).setGpuInfo(GpuConfMinimizer.Type.Cuda, 1, 2).build());
	}
	
	@Test
	public void testCudaCCDConfMinmizer1Stream() {
		check((ffparams, intergen, confSpace) -> new GpuConfMinimizer.Builder(ffparams, intergen, confSpace).setGpuInfo(GpuConfMinimizer.Type.CudaCCD, 1, 1).build());
	}
	
	@Test
	public void testCudaCCDConfMinmizer2Streams() {
		check((ffparams, intergen, confSpace) -> new GpuConfMinimizer.Builder(ffparams, intergen, confSpace).setGpuInfo(GpuConfMinimizer.Type.CudaCCD, 1, 2).build());
	}
	
	@Test
	public void testOpenCLConfMinmizer1Stream() {
		check((ffparams, intergen, confSpace) -> new GpuConfMinimizer.Builder(ffparams, intergen, confSpace).setGpuInfo(GpuConfMinimizer.Type.OpenCL, 1, 1).build());
	}
	
	@Test
	public void testOpenCLConfMinmizer2Streams() {
		check((ffparams, intergen, confSpace) -> new GpuConfMinimizer.Builder(ffparams, intergen, confSpace).setGpuInfo(GpuConfMinimizer.Type.OpenCL, 1, 2).build());
	}
	
	private static interface MinimizerFactory {
		ConfMinimizer make(ForcefieldParams ffparams, Factory<ForcefieldInteractions,Molecule> intergen, ConfSpace confSpace);
	}
	
	private void check(MinimizerFactory factory) {
		
		for (boolean doSolv : Arrays.asList(true, false)) {
			
			Info info = Infos.get(doSolv);
			ConfMinimizer minimizer = factory.make(info.ffparams, info.intergen, info.search.confSpace);
			List<EnergiedConf> econfs = minimizer.minimize(info.confs);
			minimizer.cleanup();
			
			assertThat(econfs.size(), is(info.confs.size()));
			
			for (int i=0; i<info.confs.size(); i++) {
				
				ScoredConf conf = info.confs.get(i);
				EnergiedConf econf = econfs.get(i);
				
				assertThat(econf.getAssignments(), is(conf.getAssignments()));
				assertThat(econf.getScore(), is(conf.getScore()));
				
				// penalize large errors, but not lower energies
				double absErr = econf.getEnergy() - info.expectedEnergies[i];
				assertThat(absErr, lessThanOrEqualTo(Epsilon));
			}
		}
	}
}