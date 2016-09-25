package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.math.BigDecimal;

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
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.control.MinimizingEnergyCalculator;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class TestSimplePartitionFunction extends TestBase {
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	public static SimplePartitionFunction makePfunc(SearchProblem search) {
		return makePfunc(search, 0, 0);
	}
	
	public static SimplePartitionFunction makePfunc(SearchProblem search, int numGpus, int numThreads) {
		
		// make the A* tree factory
		ConfSearchFactory confSearchFactory = new ConfSearchFactory() {
			@Override
			public ConfSearch make(EnergyMatrix emat, PruningMatrix pmat) {
				
                AStarScorer gscorer = new PairwiseGScorer(emat);
				AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), emat, 1, 0.0001);
				AStarOrder order = new StaticScoreHMeanAStarOrder();
                RCs rcs = new RCs(pmat);
                
                return new ConfAStarTree(order, gscorer, hscorer, rcs);
			}
		};
		
		// make the conf energy calculator
		ConfEnergyCalculator.Async ecalc = MinimizingEnergyCalculator.make(search, numGpus, numThreads, 0);
		
		// make the pfunc
		return new SimplePartitionFunction(search.emat, search.pruneMat, confSearchFactory, ecalc);
	}
	
	private void testProtein(int numGpus, int numThreads) {
		
		KSSearchProblem search = TestPartitionFunction.makeSearch(KSTermini.PROTEIN, "648", "654", "649 650 651 654"); 
		SimplePartitionFunction pfunc = makePfunc(search, numGpus, numThreads);

		// compute it
		final double targetEpsilon = 0.05;
		pfunc.init(targetEpsilon);
		pfunc.compute();
	
		// check the answer
		assertThat(pfunc.getStatus(), is(PartitionFunction.Status.Estimated));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		assertThat(pfunc.getValues().qstar, isRelatively(new BigDecimal("4.3704590631e+04"), targetEpsilon));
	}
	
	@Test
	public void testProtein() {
		testProtein(0, 0);
	}
	
	@Test
	public void testProteinParallel() {
		testProtein(0, 2);
	}
	
	@Test
	public void testProteinGpu() {
		testProtein(1, 0);
	}
	
	public void testLigand(int numGpus, int numThreads) {
		
		KSSearchProblem search = TestPartitionFunction.makeSearch(KSTermini.LIGAND, "155", "194", "156 172 192 193");
		SimplePartitionFunction pfunc = makePfunc(search, numGpus, numThreads);

		// compute it
		final double targetEpsilon = 0.05;
		pfunc.init(targetEpsilon);
		pfunc.compute();
	
		// check the answer
		assertThat(pfunc.getStatus(), is(PartitionFunction.Status.Estimated));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		assertThat(pfunc.getValues().qstar, isRelatively(new BigDecimal("4.4699772362e+30"), targetEpsilon));
	}
	
	@Test
	public void testLigand() {
		testLigand(0, 0);
	}
	
	@Test
	public void testLigandParallel() {
		testLigand(0, 2);
	}
	
	@Test
	public void testLigandGpu() {
		testLigand(1, 0);
	}
	
	public void testComplex(int numGpus, int numThreads) {
		
		KSSearchProblem search = TestPartitionFunction.makeSearch(KSTermini.COMPLEX, null, null, "649 650 651 654 156 172 192 193");
		SimplePartitionFunction pfunc = makePfunc(search, numGpus, numThreads);

		// compute it
		final double targetEpsilon = 0.8;
		pfunc.init(targetEpsilon);
		pfunc.compute();
	
		// check the answer
		assertThat(pfunc.getStatus(), is(PartitionFunction.Status.Estimated));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		assertThat(pfunc.getValues().qstar, isRelatively(new BigDecimal("3.5178662402e+54"), targetEpsilon));
	}
	
	@Test
	public void testComplex() {
		testComplex(0, 0);
	}
	
	@Test
	public void testComplexParallel() {
		testComplex(0, 2);
	}
	
	@Test
	public void testComplexGpu() {
		testComplex(1, 0);
	}
}
