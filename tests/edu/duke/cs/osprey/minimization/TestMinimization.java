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
	public void testCpuConfMinimizer1Thread() {
		check(new CpuConfMinimizer.Builder(ffparams, intergen, search.confSpace).build());
	}
	
	@Test
	public void testCpuConfMinimizer2Threads() {
		check(new CpuConfMinimizer.Builder(ffparams, intergen, search.confSpace).setNumThreads(2).build());
	}
	
	@Test
	public void testCudaConfMinmizer1Stream() {
		check(new GpuConfMinimizer.Builder(ffparams, intergen, search.confSpace).setGpuInfo(GpuConfMinimizer.Type.Cuda, 1, 1).build());
	}
	
	@Test
	public void testCudaConfMinmizer2Streams() {
		check(new GpuConfMinimizer.Builder(ffparams, intergen, search.confSpace).setGpuInfo(GpuConfMinimizer.Type.Cuda, 1, 2).build());
	}
	
	@Test
	public void testCudaCCDConfMinmizer1Stream() {
		check(new GpuConfMinimizer.Builder(ffparams, intergen, search.confSpace).setGpuInfo(GpuConfMinimizer.Type.CudaCCD, 1, 1).build());
	}
	
	@Test
	public void testCudaCCDConfMinmizer2Streams() {
		check(new GpuConfMinimizer.Builder(ffparams, intergen, search.confSpace).setGpuInfo(GpuConfMinimizer.Type.CudaCCD, 1, 2).build());
	}
	
	@Test
	public void testOpenCLConfMinmizer1Stream() {
		check(new GpuConfMinimizer.Builder(ffparams, intergen, search.confSpace).setGpuInfo(GpuConfMinimizer.Type.OpenCL, 1, 1).build());
	}
	
	@Test
	public void testOpenCLConfMinmizer2Streams() {
		check(new GpuConfMinimizer.Builder(ffparams, intergen, search.confSpace).setGpuInfo(GpuConfMinimizer.Type.OpenCL, 1, 2).build());
	}
	
	private void check(ConfMinimizer minimizer) {
		List<EnergiedConf> econfs = minimizer.minimize(confs);
		minimizer.cleanup();
		assertEnergies(econfs);
	}
}