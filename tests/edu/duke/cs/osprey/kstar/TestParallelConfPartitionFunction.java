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
	
	public static class Pfunc {
		
		public ConfEnergyCalculator.Async ecalc;
		public ParallelConfPartitionFunction pfunc;
		
		public Pfunc(ConfEnergyCalculator.Async ecalc, SearchProblem search, ConfSearchFactory confSearchFactory) {
			this.ecalc = ecalc;
			this.pfunc = new ParallelConfPartitionFunction(search.emat, search.pruneMat, confSearchFactory, ecalc);
		}
		
		public void cleanup() {
			this.ecalc.cleanup();
		}
	}
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	public static Pfunc makePfunc(SearchProblem search) {
		return makePfunc(search, 0, 0, 0);
	}
	
	public static Pfunc makePfunc(SearchProblem search, int numGpus, int numStreamsPerGpu, int numThreads) {
		
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
		ConfEnergyCalculator.Async ecalc = MinimizingEnergyCalculator.make(makeDefaultFFParams(), search, numGpus, numStreamsPerGpu, numThreads, true);
		
		return new Pfunc(ecalc, search, confSearchFactory);
	}
	
	private void testProtein(int numGpus, int numStreamsPerGpu, int numThreads) {
		
		KSSearchProblem search = TestPartitionFunction.makeSearch(KSTermini.PROTEIN, "648", "654", "649 650 651 654"); 
		Pfunc pfunc = makePfunc(search, numGpus, numStreamsPerGpu, numThreads);

		// compute it
		final double targetEpsilon = 0.05;
		pfunc.pfunc.init(targetEpsilon);
		pfunc.pfunc.compute();
		
		assertPfunc(pfunc.pfunc, PartitionFunction.Status.Estimated, targetEpsilon, "4.3704590631e+04" /* e=0.05 */);
		pfunc.cleanup();
	}
	
	public static void assertPfunc(PartitionFunction pfunc, PartitionFunction.Status status, double targetEpsilon, String approxQstar) {
		assertThat(pfunc.getStatus(), is(status));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		double qbound = new BigDecimal(approxQstar).doubleValue()*(1.0 - targetEpsilon);
		assertThat(pfunc.getValues().qstar.doubleValue(), greaterThanOrEqualTo(qbound));
	}
	
	@Test
	public void testProteinCpu1() {
		testProtein(0, 0, 1);
	}
	
	@Test
	public void testProteinCpu2() {
		testProtein(0, 0, 2);
	}
	
	@Test
	public void testProteinGpu1() {
		testProtein(1, 1, 0);
	}
	
	@Test
	public void testProteinGpu2() {
		testProtein(1, 2, 0);
	}
	
	public void testLigand(int numGpus, int numStreamsPerGpu, int numThreads) {
		
		KSSearchProblem search = TestPartitionFunction.makeSearch(KSTermini.LIGAND, "155", "194", "156 172 192 193");
		Pfunc pfunc = makePfunc(search, numGpus, numStreamsPerGpu, numThreads);

		// compute it
		final double targetEpsilon = 0.05;
		pfunc.pfunc.init(targetEpsilon);
		pfunc.pfunc.compute();
	
		assertPfunc(pfunc.pfunc, PartitionFunction.Status.Estimated, targetEpsilon, "4.4699772362e+30" /* e=0.05 */);
	}
	
	@Test
	public void testLigandCpu1() {
		testLigand(0, 0, 1);
	}
	
	@Test
	public void testLigandCpu2() {
		testLigand(0, 0, 2);
	}
	
	@Test
	public void testLigandGpu1() {
		testLigand(1, 1, 0);
	}
	
	@Test
	public void testLigandGpu2() {
		testLigand(1, 2, 0);
	}
	
	public void testComplex(int numGpus, int numStreamsPerGpu, int numThreads) {
		
		KSSearchProblem search = TestPartitionFunction.makeSearch(KSTermini.COMPLEX, null, null, "649 650 651 654 156 172 192 193");
		Pfunc pfunc = makePfunc(search, numGpus, numStreamsPerGpu, numThreads);

		// compute it
		final double targetEpsilon = 0.8;
		pfunc.pfunc.init(targetEpsilon);
		pfunc.pfunc.compute();
	
		assertPfunc(pfunc.pfunc, PartitionFunction.Status.Estimated, targetEpsilon, "3.5213742379e+54" /* e=0.05 */);
		pfunc.cleanup();
	}
	
	@Test
	public void testComplexCpu1() {
		testComplex(0, 0, 1);
	}
	
	@Test
	public void testComplexCpu2() {
		testComplex(0, 0, 2);
	}
	
	@Test
	public void testComplexGpu1() {
		testComplex(1, 1, 0);
	}
	
	@Test
	public void testComplexGpu2() {
		testComplex(1, 2, 0);
	}
}
