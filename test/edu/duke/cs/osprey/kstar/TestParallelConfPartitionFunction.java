package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.math.BigDecimal;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.gmec.ConfSearchFactory;
import edu.duke.cs.osprey.gmec.MinimizingConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.ParallelConfPartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class TestParallelConfPartitionFunction extends TestBase {
	
	public static class Pfunc {
		
		public MinimizingConfEnergyCalculator ecalc;
		public ParallelConfPartitionFunction pfunc;
		
		public Pfunc(MinimizingConfEnergyCalculator ecalc, SearchProblem search, ConfSearchFactory confSearchFactory) {
			this.ecalc = ecalc;
			this.pfunc = new ParallelConfPartitionFunction(search.emat, search.pruneMat, confSearchFactory, ecalc);
		}
		
		public void cleanup() {
			this.ecalc.clean();
		}
	}
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	public static Pfunc makePfunc(SearchProblem search) {
		return makePfunc(search, Parallelism.makeCpu(1));
	}
	
	public static Pfunc makePfunc(SearchProblem search, Parallelism parallelism) {
		
		// make the A* tree factory
		ConfSearchFactory confSearchFactory = new ConfSearchFactory() {
			@Override
			public ConfSearch make(EnergyMatrix emat, PruningMatrix pmat) {
				return new ConfAStarTree.Builder(emat, pmat)
					.setMPLP(new ConfAStarTree.MPLPBuilder()
						.setNumIterations(1)
					).build();
			}
		};
		
		// make the conf energy calculator
		MinimizingConfEnergyCalculator ecalc = MinimizingConfEnergyCalculator.make(makeDefaultFFParams(), search, parallelism);
		
		return new Pfunc(ecalc, search, confSearchFactory);
	}

	private void testStrand(Parallelism parallelism, KSSearchProblem search, double targetEpsilon, String approxQStar) {
		Pfunc pfunc = makePfunc(search, parallelism);
		try {
			pfunc.pfunc.init(targetEpsilon);
			pfunc.pfunc.compute();
			assertPfunc(pfunc.pfunc, PartitionFunction.Status.Estimated, targetEpsilon, approxQStar);
		} finally {
			pfunc.cleanup();
		}
	}

	public static void assertPfunc(PartitionFunction pfunc, PartitionFunction.Status status, double targetEpsilon, String approxQstar) {
		assertThat(pfunc.getStatus(), is(status));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		double qbound = new BigDecimal(approxQstar).doubleValue()*(1.0 - targetEpsilon);
		assertThat(pfunc.getValues().qstar.doubleValue(), greaterThanOrEqualTo(qbound));
	}


	private static final String PdbPath2RL0 = "examples/2RL0.kstar/2RL0.min.reduce.pdb";

	private void test2RL0Protein(Parallelism parallelism) {
		KSSearchProblem search = TestPartitionFunction.makeSearch(PdbPath2RL0, 0, "648", "654", "649 650 651 654");
		final double targetEpsilon = 0.05;
		final String approxQStar = "4.3704590631e+04"; // e=0.05
		testStrand(parallelism, search, targetEpsilon, approxQStar);
	}
	@Test public void test2RL0ProteinCpu1() { test2RL0Protein(Parallelism.makeCpu(1)); }
	@Test public void test2RL0ProteinCpu2() { test2RL0Protein(Parallelism.makeCpu(2)); }
	@Test public void test2RL0ProteinGpu1() { test2RL0Protein(Parallelism.make(4, 1, 1)); }
	@Test public void test2RL0ProteinGpu2() { test2RL0Protein(Parallelism.make(4, 1, 2)); }
	
	public void test2RL0Ligand(Parallelism parallelism) {
		KSSearchProblem search = TestPartitionFunction.makeSearch(PdbPath2RL0, 1, "155", "194", "156 172 192 193");
		final double targetEpsilon = 0.05;
		final String approxQStar = "4.4699772362e+30"; // e=0.05
		testStrand(parallelism, search, targetEpsilon, approxQStar);
	}
	@Test public void test2RL0LigandCpu1() { test2RL0Ligand(Parallelism.makeCpu(1)); }
	@Test public void test2RL0LigandCpu2() { test2RL0Ligand(Parallelism.makeCpu(2)); }
	@Test public void test2RL0LigandGpu1() { test2RL0Ligand(Parallelism.make(4, 1, 1)); }
	@Test public void test2RL0LigandGpu2() { test2RL0Ligand(Parallelism.make(4, 1, 2)); }
	
	public void test2RL0Complex(Parallelism parallelism) {
		KSSearchProblem search = TestPartitionFunction.makeSearch(PdbPath2RL0, 2, null, null, "649 650 651 654 156 172 192 193");
		final double targetEpsilon = 0.8;
		final String approxQStar = "3.5213742379e+54"; // e=0.05
		testStrand(parallelism, search, targetEpsilon, approxQStar);
	}
	@Test public void test2RL0ComplexCpu1() { test2RL0Complex(Parallelism.makeCpu(1)); }
	@Test public void test2RL0ComplexCpu2() { test2RL0Complex(Parallelism.makeCpu(2)); }
	@Test public void test2RL0ComplexGpu1() { test2RL0Complex(Parallelism.make(4, 1, 1)); }
	@Test public void test2RL0ComplexGpu2() { test2RL0Complex(Parallelism.make(4, 1, 2)); }


	private static final String PdbPath1GUA = "test-resources/1gua_adj.min.pdb";

	@Test
	public void test1GUA11Protein() {
		KSSearchProblem search = TestPartitionFunction.makeSearch(PdbPath1GUA, 0, "1", "180", "21 24 25 27 29 40");
		final double targetEpsilon = 0.1;
		final String approxQStar = "1.1838e+42"; // e=0.1
		testStrand(Parallelism.makeCpu(1), search, targetEpsilon, approxQStar);
	}

	@Test
	public void test1GUA11Ligand() {
		KSSearchProblem search = TestPartitionFunction.makeSearch(PdbPath1GUA, 1, "181", "215", "209 213");
		final double targetEpsilon = 0.1;
		final String approxQStar = "2.7098e+7"; // e=0.1
		testStrand(Parallelism.makeCpu(1), search, targetEpsilon, approxQStar);
	}

	@Test
	public void test1GUA11Complex() {
		// TODO: try null res numbers
		KSSearchProblem search = TestPartitionFunction.makeSearch(PdbPath1GUA, 1, "1", "215", "21 24 25 27 29 40 209 213");
		final double targetEpsilon = 0.1; // TEMP, set to 0.9
		final String approxQStar = "1.1195e+66"; // e=0.1
		testStrand(Parallelism.makeCpu(1), search, targetEpsilon, approxQStar);
	}
}
