package edu.duke.cs.osprey;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.HashCalculator;

public class TestAStar extends TestBase {
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	private static Map<EnergyMatrixConfig,EnergyMatrix> m_energyMatrixCache;
	
	static {
		m_energyMatrixCache = new HashMap<>();
	}
	
	private static class EnergyMatrixConfig {
		
		String pdbPath;
		int numFlexible;
		boolean addWtRots;
		boolean doMinimize;
		
		@Override
		public boolean equals(Object other) {
			if (other instanceof EnergyMatrixConfig) {
				return equals((EnergyMatrixConfig)other);
			}
			return false;
		}
		
		public boolean equals(EnergyMatrixConfig other) {
			return this.pdbPath.equals(other.pdbPath)
				&& this.numFlexible == other.numFlexible
				&& this.addWtRots == other.addWtRots
				&& this.doMinimize == other.doMinimize;
		}
		
		@Override
		public int hashCode() {
			return HashCalculator.combineHashes(
				this.pdbPath.hashCode(),
				Integer.valueOf(this.numFlexible).hashCode(),
				Boolean.valueOf(this.addWtRots).hashCode(),
				Boolean.valueOf(this.doMinimize).hashCode()
			);
		}
	}
	
	private SearchProblem makeSearchProblem(EnergyMatrixConfig emConfig) {
		
		// make the search problem
		ArrayList<String> flexRes = new ArrayList<>();
		ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
		for (int i=0; i<emConfig.numFlexible; i++) {
			flexRes.add(Integer.toString(i + 1));
			allowedAAs.add(new ArrayList<String>());
		}
		boolean addWt = true;
		boolean useEpic = false;
		boolean useTupleExpansion = false;
		boolean useEllipses = false;
		boolean useERef = false;
		boolean addResEntropy = false;
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		SearchProblem search = new SearchProblem(
			"test", emConfig.pdbPath, 
			flexRes, allowedAAs, addWt, emConfig.doMinimize, useEpic, new EPICSettings(), useTupleExpansion,
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, emConfig.addWtRots
		);
		
		// calculate the energy matrix, but check the cache first
		search.emat = m_energyMatrixCache.get(emConfig);
		if (search.emat == null) {
			EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(search.confSpace, search.shellResidues, useERef, addResEntropy);
			emCalc.calcPEM();
			search.emat = emCalc.getEMatrix();
			m_energyMatrixCache.put(emConfig, search.emat);
		}
		
		// calculate an "identity" pruning matrix (ie, no pruning)
		search.pruneMat = new PruningMatrix(search.confSpace, search.emat.getPruningInterval());
		
		return search;
	}
	
	
	// RIGID TESTS (without pruning)
	
	private SearchProblem makeSearchProblemDagkRigid() {
		EnergyMatrixConfig emConfig = new EnergyMatrixConfig();
		emConfig.pdbPath = "test/DAGK/2KDC.P.forOsprey.pdb";
		emConfig.numFlexible = 8;
		emConfig.addWtRots = true;
		emConfig.doMinimize = false;
		return makeSearchProblem(emConfig);
	}
	
	private void checkDagkRigid(ConfSearch tree, SearchProblem search) {
		int[] conf = tree.nextConf();
		double confEnergy = search.emat.getInternalEnergy(new RCTuple(conf));
	
		assertThat(confEnergy, isRelatively(-78.78903548260008));
		assertThat(conf, is(new int[] { 0, 6, 7, 0, 16, 1, 0, 6 }));
	}
	
	@Test
	public void testDagkRigidConfTree() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		ConfSearch tree = ConfTree.makeFull(search);
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidDynamicScoreOrderTraditional() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new DynamicHMeanAStarOrder();
		AStarScorer hscorer = new TraditionalPairwiseHScorer(search.emat, rcs);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidStaticScoreOrderTraditional() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new TraditionalPairwiseHScorer(search.emat, rcs);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidDynamicOrderMPLP0Iter() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new DynamicHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(null, search.emat, 0, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidStaticScoreOrderMPLP0Iter() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(null, search.emat, 0, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidStaticScoreOrderMPLPNode1Iter() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidStaticScoreOrderMPLPNode5Iter() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 5, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testRigidStaticScoreOrderMPLPEdge1Iter() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new EdgeUpdater(), search.emat, 1, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testRigidStaticScoreOrderMPLPEdge20Iter() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new EdgeUpdater(), search.emat, 20, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidDynamicOrderMPLPNode1Iter() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new DynamicHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkRigid(tree, search);
	}
	
	
	// RIGID TESTS (with pruning)
	
	private void setPruning(SearchProblem search, int[] partialConf) {
		
		int numChosen = 0;
		for (int pos=0; pos<search.confSpace.numPos; pos++) {
			int chosenRc = partialConf[pos];
			if (chosenRc > 0) {
				numChosen++;
				
				// prune all but this rc
				for (int rc=0; rc<search.confSpace.posFlex.get(pos).RCs.size(); rc++) {
					if (rc != chosenRc) {
						search.pruneMat.setOneBody(pos, rc, true);
					}
				}
			}
		}
		
		assert (new RCs(search.pruneMat).getNumTrivialPos() == numChosen);
	}
	
	private SearchProblem makePrunedSearchProblemDagkRigid() {
		SearchProblem search = makeSearchProblemDagkRigid();
		setPruning(search, new int[] { -1, -1, 7, -1, 16, -1, 0, 6 });
		return search;
	}
	
	@Test
	public void testDagkRigidConfTreePruned() {
		SearchProblem search = makePrunedSearchProblemDagkRigid();
		
		ConfSearch tree = ConfTree.makeFull(search);
		
		checkDagkRigid(tree, search);
	}

	@Test
	public void testDagkRigidStaticScoreOrderMPLPNode1IterPruned() {
		SearchProblem search = makePrunedSearchProblemDagkRigid();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkRigid(tree, search);
	}
	
	
	// CONTINUOUS TESTS
	
	private SearchProblem makeSearchProblemDagkContinuous() {
		EnergyMatrixConfig emConfig = new EnergyMatrixConfig();
		emConfig.pdbPath = "test/DAGK/2KDC.P.forOsprey.pdb";
		emConfig.numFlexible = 8;
		emConfig.addWtRots = true;
		emConfig.doMinimize = true;
		return makeSearchProblem(emConfig);
	}
	
	private void checkDagkContinuous(ConfSearch tree, SearchProblem search) {
		int[] conf = tree.nextConf();
		double confEnergy = search.emat.getInternalEnergy(new RCTuple(conf));
	
		assertThat(confEnergy, isRelatively(-84.85105599013883));
		assertThat(conf, is(new int[] { 0, 0, 7, 16, 16, 1, 1, 17 }));
	}
	
	@Test
	public void testDagkContinuousConfTree() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		ConfSearch tree = ConfTree.makeFull(search);
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousDynamicScoreOrderTraditional() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new DynamicHMeanAStarOrder();
		AStarScorer hscorer = new TraditionalPairwiseHScorer(search.emat, rcs);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousStaticScoreOrderTraditional() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new TraditionalPairwiseHScorer(search.emat, rcs);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousDynamicOrderMPLP0Iter() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new DynamicHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(null, search.emat, 0, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousStaticScoreOrderMPLP0Iter() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(null, search.emat, 0, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousStaticScoreOrderMPLPNode1Iter() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousStaticScoreOrderMPLPNode5Iter() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 5, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testContinuousStaticScoreOrderMPLPEdge1Iter() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new EdgeUpdater(), search.emat, 1, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testContinuousStaticScoreOrderMPLPEdge20Iter() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new EdgeUpdater(), search.emat, 20, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousDynamicOrderMPLPNode1Iter() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new DynamicHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		
		checkDagkContinuous(tree, search);
	}
}
