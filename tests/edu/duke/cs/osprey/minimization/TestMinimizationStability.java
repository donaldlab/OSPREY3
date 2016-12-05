package edu.duke.cs.osprey.minimization;

import static org.junit.Assert.*;

import java.io.File;
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
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class TestMinimizationStability extends TestBase {
	
	private static final double EnergyEpsilon = 1e-8;
	
	private static SearchProblem search;
	private static List<ScoredConf> confs;
	private static EnergyFunctionGenerator egen;
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
		
		// for these small problems, more than one thread is actually slower
		MultiTermEnergyFunction.setNumThreads(1);
		
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
		
		search = new SearchProblem(
			"test", "test/1CC8/1CC8.ss.pdb", 
			resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null,
			false, new ArrayList<>()
		);
		
		egen = new EnergyFunctionGenerator(makeDefaultFFParams(), Double.POSITIVE_INFINITY, false);
		
		// calc the energy matrix
		File ematFile = new File("/tmp/testMinimizationStability.emat.dat");
		if (ematFile.exists()) {
			search.emat = (EnergyMatrix)ObjectIO.readObject(ematFile.getAbsolutePath(), false);
		} else {
			ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
			tasks.start(2);
			SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, search.confSpace, search.shellResidues);
			search.emat = new SimpleEnergyMatrixCalculator(ecalc).calcEnergyMatrix(tasks);
			tasks.stop();
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
		}
		
		// don't prune anything
		search.pruneMat = new PruningMatrix(search.confSpace, 1000);
		
		// build A* tree
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		// get the confs
		final int NumConfs = 10;
		confs = new ArrayList<>();
		for (int i=0; i<NumConfs; i++) {
			confs.add(tree.nextConf());
		}
	}
	
	@Test
	public void sharedMol() {
		
		ParameterizedMoleculeCopy pmol = new ParameterizedMoleculeCopy(search.confSpace);
		
		System.out.println("minimizing in ascending order...");
		double[] ascendingEnergies = new double[confs.size()];
		for (int i=0; i<confs.size(); i++) {
			ascendingEnergies[i] = minimize(confs.get(i), pmol);
		}
		
		System.out.println("minimizing in descending order...");
		double[] descendingEnergies = new double[confs.size()];
		for (int i=confs.size() - 1; i>=0; i--) {
			descendingEnergies[i] = minimize(confs.get(i), pmol);
		}
		
		checkResults(ascendingEnergies, descendingEnergies, EnergyEpsilon);
	}
	
	@Test
	public void separateMols() {
		
		System.out.println("minimizing in ascending order...");
		double[] ascendingEnergies = new double[confs.size()];
		for (int i=0; i<confs.size(); i++) {
			ParameterizedMoleculeCopy pmol = new ParameterizedMoleculeCopy(search.confSpace);
			ascendingEnergies[i] = minimize(confs.get(i), pmol);
		}
		
		System.out.println("minimizing in descending order...");
		double[] descendingEnergies = new double[confs.size()];
		for (int i=confs.size() - 1; i>=0; i--) {
			ParameterizedMoleculeCopy pmol = new ParameterizedMoleculeCopy(search.confSpace);
			descendingEnergies[i] = minimize(confs.get(i), pmol);
		}
		
		// energies should be exactly the same
		checkResults(ascendingEnergies, descendingEnergies, 0);
	}
	
	private double minimize(ScoredConf conf, ParameterizedMoleculeCopy pmol) {
		EnergyFunction efunc = egen.fullConfEnergy(search.confSpace, search.shellResidues, pmol.getCopiedMolecule());
		return new ConfMinimizer().minimize(pmol, conf, efunc, search.confSpace).getEnergy();
	}
	
	private void checkResults(double[] ascendingEnergies, double[] descendingEnergies, double energyEpsilon) {
		
		double sum = 0;
		
		int n = ascendingEnergies.length;
		for (int i=0; i<n; i++) {
			
			double absErr = Math.abs(ascendingEnergies[i] - descendingEnergies[i]);
			sum += absErr*absErr;
			
			System.out.println(String.format("asc: %20.12f   desc: %20.12f   err: %e   %s",
				ascendingEnergies[i], descendingEnergies[i],
				absErr,
				absErr > EnergyEpsilon ? "   <-- too much err" : ""
			));
		}
		
		double rmsd = Math.sqrt(sum/n);
		System.out.println(String.format("RMS err: %12.12f", rmsd));
		
		// check the results
		for (int i=0; i<n; i++) {
			assertThat(ascendingEnergies[i], isAbsolutely(descendingEnergies[i], energyEpsilon));
		}
	}
	
}
