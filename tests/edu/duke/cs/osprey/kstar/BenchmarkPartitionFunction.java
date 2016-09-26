package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.math.BigDecimal;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;
import edu.duke.cs.osprey.tools.Stopwatch;

public class BenchmarkPartitionFunction extends TestBase {
	
	private static final int NumThreads = 2;
	
	public static void main(String[] args) {
		
		initDefaultEnvironment();
		
		// disable energy function-level parallelism
		// it's actually slower for small designs
		MultiTermEnergyFunction.setNumThreads(1);
		
		// set fork join pool parallelism used by PFParallelN
		ThreadParallelism.setNumThreads(NumThreads);
		
		benchmarkProtein();
		benchmarkLigand();
		benchmarkComplex();
	}
	
	private static void benchmark(KSSearchProblem search, int strand, double targetEpsilon, String qstar) {
		
		final boolean reportProgress = false;
		
		// setup the config for the pffactory
		KSConfigFileParser cfp = new KSConfigFileParser();
		cfp.getParams().setValue("MinimizationThreads", Integer.toString(NumThreads));
		
		System.out.println("\n\nBenchmarking " + KSTermini.getTerminiString(strand) + "...\n");
		
		// test parallel implementation
		PFAbstract pfunc = TestPartitionFunction.makePfunc(search, "parallel0", KSTermini.PROTEIN, null, cfp);
		PFAbstract.suppressOutput = !reportProgress;
		PFAbstract.targetEpsilon = targetEpsilon;
		
		System.out.println("computing pfunc " + pfunc.getClass().getSimpleName() + " ...");
		Stopwatch stopwatchParallel = new Stopwatch().start();
		pfunc.start();
		pfunc.runToCompletion();
		System.out.println(String.format("finished in %s", stopwatchParallel.stop().getTime(2)));
		
		// test simple implementation
		SimplePartitionFunction spfunc = TestSimplePartitionFunction.makePfunc(search, 0, NumThreads);
		spfunc.init(targetEpsilon);
		spfunc.setReportProgress(reportProgress);
		
		System.out.println("computing pfunc " + spfunc.getClass().getSimpleName() + " ...");
		Stopwatch stopwatchSimple = new Stopwatch().start();
		spfunc.compute();
		System.out.println(String.format("finished in %s, speedup=%.2f", stopwatchSimple.stop().getTime(2), (double)stopwatchParallel.getTimeNs()/stopwatchSimple.getTimeNs()));
		
		// test simple implementation on gpu
		SimplePartitionFunction spfuncgpu = TestSimplePartitionFunction.makePfunc(search, 1, NumThreads);
		spfuncgpu.init(targetEpsilon);
		spfuncgpu.setReportProgress(reportProgress);
		
		System.out.println("computing pfunc " + spfuncgpu.getClass().getSimpleName() + " on GPU ...");
		Stopwatch stopwatchSimpleGpu = new Stopwatch().start();
		spfuncgpu.compute();
		System.out.println(String.format("finished in %s, speedup=%.2f", stopwatchSimpleGpu.stop().getTime(2), (double)stopwatchParallel.getTimeNs()/stopwatchSimpleGpu.getTimeNs()));
		
		// test adapted simple implementation
		PFAbstract pfuncAdapted = TestPartitionFunction.makePfunc(search, "simple", KSTermini.PROTEIN, null, cfp);
		
		System.out.println("computing pfunc " + pfuncAdapted.getClass().getSimpleName() + " ...");
		Stopwatch stopwatchAdapted = new Stopwatch().start();
		pfuncAdapted.start();
		pfuncAdapted.runToCompletion();
		System.out.println(String.format("finished in %s, speedup=%.2f", stopwatchAdapted.stop().getTime(2), (double)stopwatchParallel.getTimeNs()/stopwatchAdapted.getTimeNs()));
		
		// check the results, just in case
		checkProteinPfunc(pfunc, targetEpsilon, qstar);
		checkProteinPfunc(spfunc, targetEpsilon, qstar);
		checkProteinPfunc(spfuncgpu, targetEpsilon, qstar);
		checkProteinPfunc(pfuncAdapted, targetEpsilon, qstar);
		
		System.out.println();
	}
	
	private static void checkProteinPfunc(PFAbstract pfunc, double targetEpsilon, String qstar) {
		assertThat(pfunc.getEpsilonStatus(), is(EApproxReached.TRUE));
		assertThat(pfunc.getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		assertThat(pfunc.getQStar(), isRelatively(new BigDecimal(qstar), targetEpsilon));
	}
	
	private static void checkProteinPfunc(PartitionFunction pfunc, double targetEpsilon, String qstar) {
		assertThat(pfunc.getStatus(), is(PartitionFunction.Status.Estimated));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		assertThat(pfunc.getValues().qstar, isRelatively(new BigDecimal(qstar), targetEpsilon));
	}
	
	private static void benchmarkProtein() {
		
		final double targetEpsilon = 0.01;
		final String qstar = "4.3704590631e+04";
		int strand = KSTermini.PROTEIN;
		
		KSSearchProblem search = TestPartitionFunction.makeSearch(strand, "648", "654", "649 650 651 654");
		
		benchmark(search, strand, targetEpsilon, qstar);
	}
	
	private static void benchmarkLigand() {
		
		final double targetEpsilon = 0.01;
		final String qstar = "4.4699772362e+30";
		int strand = KSTermini.LIGAND;
		
		KSSearchProblem search = TestPartitionFunction.makeSearch(strand, "155", "194", "156 172 192 193");
		
		benchmark(search, strand, targetEpsilon, qstar);
	}
	
	private static void benchmarkComplex() {
		
		final double targetEpsilon = 0.8;
		final String qstar = "3.5178662402e+54"; 
		int strand = KSTermini.COMPLEX;
		
		KSSearchProblem search = TestPartitionFunction.makeSearch(strand, null, null, "649 650 651 654 156 172 192 193");
		
		benchmark(search, strand, targetEpsilon, qstar);
	}
}
