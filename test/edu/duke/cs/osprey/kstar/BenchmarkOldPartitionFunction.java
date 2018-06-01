/*
 ** This file is part of OSPREY 3.0
 **
 ** OSPREY Protein Redesign Software Version 3.0
 ** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
 **
 ** OSPREY is free software: you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License version 2
 ** as published by the Free Software Foundation.
 **
 ** You should have received a copy of the GNU General Public License
 ** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
 **
 ** OSPREY relies on grants for its development, and since visibility
 ** in the scientific literature is essential for our success, we
 ** ask that users of OSPREY cite our papers. See the CITING_OSPREY
 ** document in this distribution for more information.
 **
 ** Contact Info:
 **    Bruce Donald
 **    Duke University
 **    Department of Computer Science
 **    Levine Science Research Center (LSRC)
 **    Durham
 **    NC 27708-0129
 **    USA
 **    e-mail: www.cs.duke.edu/brd/
 **
 ** <signature of Bruce Donald>, Mar 1, 2018
 ** Bruce Donald, Professor of Computer Science
 */

package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.math.BigDecimal;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.TestParallelConfPartitionFunction.Pfunc;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;
import edu.duke.cs.osprey.tools.Stopwatch;

public class BenchmarkOldPartitionFunction extends TestBase {

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

	private static void benchmark(KSSearchProblem search, int strand, String flexibility, double targetEpsilon, String qstar) {

		final boolean reportProgress = false;

		// setup the config for the pffactory
		KSConfigFileParser cfp = new KSConfigFileParser();
		cfp.params.setValue("MinimizationThreads", Integer.toString(NumThreads));

		System.out.println("\n\nBenchmarking " + "Strand"+strand + "...\n");

		// test parallel implementation
		PFAbstract pfunc = TestPartitionFunction.makePfunc(search, "parallel0", 0, flexibility, cfp);
		PFAbstract.suppressOutput = !reportProgress;
		PFAbstract.targetEpsilon = targetEpsilon;

		System.out.println("computing pfunc " + pfunc.getClass().getSimpleName() + " ...");
		Stopwatch stopwatchParallel = new Stopwatch().start();
		pfunc.start();
		pfunc.runToCompletion();
		System.out.println(String.format("finished in %s", stopwatchParallel.stop().getTime(2)));

		// test parallel conf implementation
		Pfunc pcpfunc = TestParallelConfPartitionFunction.makePfunc(search, Parallelism.makeCpu(NumThreads));
		pcpfunc.pfunc.init(null, null, targetEpsilon);
		pcpfunc.pfunc.setReportProgress(reportProgress);

		System.out.println("computing pfunc " + pcpfunc.getClass().getSimpleName() + " ...");
		Stopwatch stopwatchSimple = new Stopwatch().start();
		pcpfunc.pfunc.compute();
		System.out.println(String.format("finished in %s, speedup=%.2f", stopwatchSimple.stop().getTime(2), (double)stopwatchParallel.getTimeNs()/stopwatchSimple.getTimeNs()));
		pcpfunc.cleanup();

		// test parallel conf implementation on gpu
		Pfunc pcpfuncgpu = TestParallelConfPartitionFunction.makePfunc(search, Parallelism.make(4, 1, 1));
		pcpfuncgpu.pfunc.init(null, null, targetEpsilon);
		pcpfuncgpu.pfunc.setReportProgress(reportProgress);

		System.out.println("computing pfunc " + pcpfuncgpu.getClass().getSimpleName() + " on GPU with 1 stream ...");
		Stopwatch stopwatchSimpleGpu = new Stopwatch().start();
		pcpfuncgpu.pfunc.compute();
		System.out.println(String.format("finished in %s, speedup=%.2f", stopwatchSimpleGpu.stop().getTime(2), (double)stopwatchParallel.getTimeNs()/stopwatchSimpleGpu.getTimeNs()));
		pcpfuncgpu.cleanup();

		// test parallel conf implementation on gpu
		Pfunc pcpfuncgpuMulti = TestParallelConfPartitionFunction.makePfunc(search, Parallelism.make(4, 1, 16));
		pcpfuncgpuMulti.pfunc.init(null, null, targetEpsilon);
		pcpfuncgpuMulti.pfunc.setReportProgress(reportProgress);

		System.out.println("computing pfunc " + pcpfuncgpuMulti.getClass().getSimpleName() + " on GPU with 16 streams ...");
		Stopwatch stopwatchSimpleGpuMulti = new Stopwatch().start();
		pcpfuncgpuMulti.pfunc.compute();
		System.out.println(String.format("finished in %s, speedup=%.2f", stopwatchSimpleGpuMulti.stop().getTime(2), (double)stopwatchParallel.getTimeNs()/stopwatchSimpleGpuMulti.getTimeNs()));
		pcpfuncgpuMulti.cleanup();

		// test adapted parallel conf implementation
		PFAbstract pfuncAdapted = TestPartitionFunction.makePfunc(search, "parallelConf", 0, flexibility, cfp);

		System.out.println("computing pfunc " + pfuncAdapted.getClass().getSimpleName() + " ...");
		Stopwatch stopwatchAdapted = new Stopwatch().start();
		pfuncAdapted.start();
		pfuncAdapted.runToCompletion();
		System.out.println(String.format("finished in %s, speedup=%.2f", stopwatchAdapted.stop().getTime(2), (double)stopwatchParallel.getTimeNs()/stopwatchAdapted.getTimeNs()));

		// check the results, just in case
		checkProteinPfunc(pfunc, targetEpsilon, qstar);
		checkProteinPfunc(pcpfunc.pfunc, targetEpsilon, qstar);
		checkProteinPfunc(pcpfuncgpu.pfunc, targetEpsilon, qstar);
		checkProteinPfunc(pcpfuncgpuMulti.pfunc, targetEpsilon, qstar);
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

	private static final String PdbPath2RL0 = "examples/2RL0.kstar/2RL0.min.reduce.pdb";

	private static void benchmarkProtein() {

		final double targetEpsilon = 0.01;
		final String qstar = "4.3704590631e+04";
		int strand = 0;
		String flexibility = "649 650 651 654";

		KSSearchProblem search = TestPartitionFunction.makeSearch(PdbPath2RL0, strand, "648", "654", flexibility);

		benchmark(search, strand, flexibility, targetEpsilon, qstar);
	}

	private static void benchmarkLigand() {

		final double targetEpsilon = 0.01;
		final String qstar = "4.4699772362e+30";
		int strand = 1;
		String flexibility = "156 172 192 193";

		KSSearchProblem search = TestPartitionFunction.makeSearch(PdbPath2RL0, strand, "155", "194", flexibility);

		benchmark(search, strand, flexibility, targetEpsilon, qstar);
	}

	private static void benchmarkComplex() {

		final double targetEpsilon = 0.1;
		final String qstar = "3.5178662402e+54";
		int strand = 2;
		String flexibility = "649 650 651 654 156 172 192 193";

		KSSearchProblem search = TestPartitionFunction.makeSearch(PdbPath2RL0, strand, null, null, flexibility);

		benchmark(search, strand, flexibility, targetEpsilon, qstar);
	}
}
