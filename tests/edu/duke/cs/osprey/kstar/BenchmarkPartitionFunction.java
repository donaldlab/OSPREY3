package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.math.BigDecimal;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;
import edu.duke.cs.osprey.tools.Stopwatch;

public class BenchmarkPartitionFunction extends TestBase {
	
	public static void main(String[] args) {
		
		initDefaultEnvironment();
		
		// disable energy function-level parallelism
		MultiTermEnergyFunction.setNumThreads(1);
		
		// cut down on spam for benchmarking
		PFAbstract.suppressOutput = true;
		
		benchmarkProtein();
	}
	
	private static void benchmarkProtein() {
		
		final double targetEpsilon = 0.01;
		
		KSSearchProblem search = TestPartitionFunction.makeSearch(KSTermini.PROTEIN, "648", "654", "649 650 651 654");
		
		System.out.println("\n");
		
		// test parallel implementation
		PFAbstract pfunc = TestPartitionFunction.makePfunc(search, "parallel0", KSTermini.PROTEIN);
		ThreadParallelism.setNumThreads(1);
		PFAbstract.targetEpsilon = targetEpsilon;
		
		Stopwatch stopwatchParallel = new Stopwatch().start();
		pfunc.start();
		pfunc.runToCompletion();
		System.out.println("Pfunc finished in " + stopwatchParallel.stop().getTime(2));
		
		checkProteinPfunc(pfunc);
		
		// test simple implementation
		SimplePartitionFunction spfunc = TestSimplePartitionFunction.makePfunc(search);
		spfunc.init(targetEpsilon);
		
		Stopwatch stopwatchSimple = new Stopwatch().start();
		spfunc.compute();
		System.out.println("Pfunc finished in " + stopwatchSimple.stop().getTime(2));
		System.out.println(String.format("\tSpeedup: %.2fx", (double)stopwatchParallel.getTimeNs()/stopwatchSimple.getTimeNs()));
		
		checkProteinPfunc(pfunc);
	}
	
	private static void checkProteinPfunc(PFAbstract pfunc) {
		assertThat(pfunc.getEpsilonStatus(), is(EApproxReached.TRUE));
		assertThat(pfunc.getEffectiveEpsilon(), lessThanOrEqualTo(PFAbstract.targetEpsilon));
		assertThat(pfunc.getQStar(), isRelatively(new BigDecimal("4.3704590631e+04"), PFAbstract.targetEpsilon));
	}
}
