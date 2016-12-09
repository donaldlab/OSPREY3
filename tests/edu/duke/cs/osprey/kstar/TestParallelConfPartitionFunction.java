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
import edu.duke.cs.osprey.kstar.pfunc.ParallelConfPartitionFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class TestParallelConfPartitionFunction extends TestBase {
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	public static ParallelConfPartitionFunction makePfunc(SearchProblem search) {
		return makePfunc(search, 0, 0);
	}
	
	public static ParallelConfPartitionFunction makePfunc(SearchProblem search, int numGpus, int numThreads) {
		
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
		ConfEnergyCalculator.Async ecalc = MinimizingEnergyCalculator.make(makeDefaultFFParams(), search, numGpus, numThreads, 0, true);
		
		// make the pfunc
		return new ParallelConfPartitionFunction(search.emat, search.pruneMat, confSearchFactory, ecalc);
	}
	
	private void testProtein(int numGpus, int numThreads) {
		
		KSSearchProblem search = TestPartitionFunction.makeSearch(KSTermini.PROTEIN, "648", "654", "649 650 651 654"); 
		ParallelConfPartitionFunction pfunc = makePfunc(search, numGpus, numThreads);

		// compute it
		final double targetEpsilon = 0.05;
		pfunc.init(targetEpsilon);
		pfunc.compute();
		
		assertPfunc(pfunc, PartitionFunction.Status.Estimated, targetEpsilon, "4.3704590631e+04" /* e=0.05 */);
	}
	
	public static void assertPfunc(PartitionFunction pfunc, PartitionFunction.Status status, double targetEpsilon, String approxQstar) {
		assertThat(pfunc.getStatus(), is(status));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		double qbound = new BigDecimal(approxQstar).doubleValue()*(1.0 - targetEpsilon);
		assertThat(pfunc.getValues().qstar.doubleValue(), greaterThanOrEqualTo(qbound));
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
		ParallelConfPartitionFunction pfunc = makePfunc(search, numGpus, numThreads);

		// compute it
		final double targetEpsilon = 0.05;
		pfunc.init(targetEpsilon);
		pfunc.compute();
	
		assertPfunc(pfunc, PartitionFunction.Status.Estimated, targetEpsilon, "4.4699772362e+30" /* e=0.05 */);
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
		ParallelConfPartitionFunction pfunc = makePfunc(search, numGpus, numThreads);

		// compute it
		final double targetEpsilon = 0.8;
		pfunc.init(targetEpsilon);
		pfunc.compute();
	
		assertPfunc(pfunc, PartitionFunction.Status.Estimated, targetEpsilon, "3.5213742379e+54" /* e=0.05 */);
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
