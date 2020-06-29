package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.coffee.directors.PfuncDirector;
import edu.duke.cs.osprey.confspace.ConfSpaceIteration;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.confspace.compiled.TestConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.compiled.EmatCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculatorAdapter;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.MathTools.*;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.slf4j.bridge.SLF4JBridgeHandler;

import java.util.function.Supplier;

import static edu.duke.cs.osprey.tools.Log.log;


public class BenchmarkCoffee {

	public static void main(String[] args) {

		// configure hazelcast logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();

		// load a complex state
		var complexClassic = TestConfSpace.Design2RL0Interface7Mut.makeClassic().complex;
		var complexCompiled = TestConfSpace.Design2RL0Interface7Mut.makeCompiled().complex;

		var seqClassic = complexClassic.seqSpace.makeWildTypeSequence();
		var seqCompiled = complexCompiled.seqSpace.makeWildTypeSequence();

		var oneCpu = Parallelism.makeCpu(1);
		var allCpus = Parallelism.makeCpu(Parallelism.getMaxNumCPUs());

		var bounds = Bounds.Classic;
		// var bonuds = Bounds.Tighter;

		var noStaticStatic = false;
		var yesStaticStatic = true;

		double epsilon = 0.68; // amazingly, also corresponds to roughly a 0.67 width in free energy
		double gWidthMax = 0.67;

		// benchmarks, on jerry4 (up to 48 threads, 24 cores)
		var cpus = Parallelism.makeCpu(6);

		//benchmark("GD     classic   noSS", () -> gradientDescent(complexClassic, seqClassic, cpus, bounds, noStaticStatic, epsilon));

		// 1 thread
		//[28    4     0     28    13    28    10   ] scores:   16421, confs: 262, score: -126.800794, energy: -120.704369, bounds:[   94.681678,   95.174941] (log10p1), delta:0.678828, time:    1.32 m, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     classic   noSS   emat   19095 ms ( 19.10 s)   pfunc   79191 ms (  1.32 m)   G [-129.9660,-129.2924]  w =  0.6736

		//benchmark("GD     compiled yesSS", () -> gradientDescent(complexCompiled, seqCompiled, cpus, bounds, yesStaticStatic, epsilon));

		// 1 thread
		//[14    40    5     42    28    13    13   ] scores:   16658, confs: 167, score:-1557.838995, energy:-1553.935900, bounds:[ 1142.602893, 1143.096744] (log10p1), delta:0.679263, time:   38.90 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     compiled yesSS   emat   11018 ms ( 11.02 s)   pfunc   38952 ms ( 38.95 s)   G [-1560.9540,-1560.2796]  w =  0.6743

		// 2 threads
		//[28    38    28    39    28    13    9    ] scores:   13754, confs: 169, score:-1557.822975, energy:-1555.779319, bounds:[ 1142.603402, 1143.092819] (log10p1), delta:0.675972, time:   21.11 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     compiled yesSS   emat   10991 ms ( 10.99 s)   pfunc   21030 ms ( 21.03 s)   G [-1560.9486,-1560.2803]  w =  0.6683

		// 6 threads
		//[14    38    28    41    28    13    9    ] scores:   13659, confs: 174, score:-1557.808425, energy:-1555.769220, bounds:[ 1142.604120, 1143.082060] (log10p1), delta:0.667295, time:    7.60 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     compiled yesSS   emat   10883 ms ( 10.88 s)   pfunc    7530 ms (  7.53 s)   G [-1560.9340,-1560.2813]  w =  0.6527

		// 12 threads
		//[14    38    28    41    28    13    9    ] scores:   13338, confs: 180, score:-1557.808425, energy:-1555.769220, bounds:[ 1142.604855, 1143.069136] (log10p1), delta:0.656664, time:    4.50 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     compiled yesSS   emat   11087 ms ( 11.09 s)   pfunc    4427 ms (  4.43 s)   G [-1560.9163,-1560.2823]  w =  0.6340

		// 24 threads
		//[14    36    28    41    28    13    13   ] scores:   16507, confs: 190, score:-1557.750911, energy:-1555.812669, bounds:[ 1142.606579, 1143.047615] (log10p1), delta:0.637788, time:    2.97 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     compiled yesSS   emat   11235 ms ( 11.24 s)   pfunc    2891 ms (  2.89 s)   G [-1560.8869,-1560.2847]  w =  0.6023

		// 48 threads
		//[14    38    28    39    28    13    9    ] scores:   10956, confs: 216, score:-1557.696283, energy:-1555.650422, bounds:[ 1142.609931, 1142.998349] (log10p1), delta:0.591133, time:    2.98 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     compiled yesSS   emat   11022 ms ( 11.02 s)   pfunc    2901 ms (  2.90 s)   G [-1560.8196,-1560.2892]  w =  0.5304

		benchmark("COFFEE compiled yesSS", () -> coffee(complexCompiled, seqCompiled, cpus, bounds, yesStaticStatic, gWidthMax));

		// 1 thread
		//COFFEE-0: 	G [-1560.937,-1560.278]   width 0.659032 of 0.000000   confs       254   avgap 2.84   nodedb  26.6%   rr Infinity   time 59.07 s
		//COFFEE compiled yesSS   emat   24947 ms ( 24.95 s)   pfunc   59074 ms ( 59.07 s)   G [-1560.9373,-1560.2783]  w =  0.6590

		// 2 threads
		//COFFEE-0: 	G [-1560.935,-1560.276]   width 0.658562 of 0.000000   confs       258   avgap 2.82   nodedb  24.2%   rr Infinity   time 33.01 s
		//COFFEE compiled yesSS   emat   11634 ms ( 11.63 s)   pfunc   33010 ms ( 33.01 s)   G [-1560.9346,-1560.2761]  w =  0.6586

		// 6 threads
		//COFFEE-0: 	G [-1560.821,-1560.281]   width 0.539462 of 0.000000   confs       273   avgap 2.73   nodedb  22.7%   rr Infinity   time 14.06 s
		//COFFEE compiled yesSS   emat    7730 ms (  7.73 s)   pfunc   14058 ms ( 14.06 s)   G [-1560.8206,-1560.2811]  w =  0.5395

		// 12 threads
		//COFFEE-0: 	G [-1560.795,-1560.287]   width 0.507800 of 0.000000   confs       277   avgap 2.72   nodedb  21.1%   rr Infinity   time 7.84 s
		//COFFEE compiled yesSS   emat    7469 ms (  7.47 s)   pfunc    7842 ms (  7.84 s)   G [-1560.7945,-1560.2867]  w =  0.5078

		// 24 threads
		//COFFEE-0: 	G [-1560.858,-1560.290]   width 0.567334 of 0.000000   confs       226   avgap 2.64   nodedb  18.8%   rr Infinity   time 4.15 s
		//COFFEE compiled yesSS   emat    7494 ms (  7.49 s)   pfunc    4147 ms (  4.15 s)   G [-1560.8576,-1560.2903]  w =  0.5673

		// 48 threads
		//COFFEE-0: 	G [-1560.680,-1560.297]   width 0.383059 of 0.000000   confs       297   avgap 2.73   nodedb  21.9%   rr Infinity   time 4.23 s
		//COFFEE compiled yesSS   emat    8059 ms (  8.06 s)   pfunc    4228 ms (  4.23 s)   G [-1560.6802,-1560.2972]  w =  0.3831

		// TODO: NEXTTIME: scale up, GPUs
		// TODO: NEXTTIME: scale out
	}

	private enum Bounds {

		Classic(EnergyPartition.Traditional, PosInterDist.DesmetEtAl1992),
		Tighter(EnergyPartition.AllOnPairs, PosInterDist.TighterBounds);

		public final EnergyPartition epart;
		public final PosInterDist posInterDist;

		Bounds(EnergyPartition epart, PosInterDist posInterDist) {
			this.epart = epart;
			this.posInterDist = posInterDist;
		}
	}

	private static class Result {
		DoubleBounds freeEnergy;
		Stopwatch pfuncStopwatch;
		Stopwatch ematStopwatch;
	}

	private static void benchmark(String name, Supplier<Result> task) {

		// run the task
		var result = task.get();

		log("%20s   emat  %6d ms (%8s)   pfunc  %6d ms (%8s)   G [%9.4f,%9.4f]  w = %7.4f",
			name,
			(int)result.ematStopwatch.getTimeMs(), result.ematStopwatch.getTime(2),
			(int)result.pfuncStopwatch.getTimeMs(), result.pfuncStopwatch.getTime(2),
			result.freeEnergy.lower, result.freeEnergy.upper, result.freeEnergy.size()
		);
	}

	private static Result gradientDescent(ConfSpaceIteration confSpace, Sequence seq, Parallelism parallelism, Bounds bounds, boolean includeStaticStatic, double epsilon) {
		var result = new Result();
		if (confSpace instanceof SimpleConfSpace) {
			if (includeStaticStatic) {
				throw new IllegalArgumentException("classic doesn't support static-static energies");
			}
			gradientDescentClassic((SimpleConfSpace)confSpace, seq, parallelism, bounds.epart, epsilon, result);
		} else if (confSpace instanceof ConfSpace) {
			gradientDescentCompiled((ConfSpace)confSpace, seq, parallelism, bounds.posInterDist, includeStaticStatic, epsilon, result);
		} else {
			throw new Error(confSpace.getClass().getName());
		}
		return result;
	}

	private static void gradientDescentClassic(SimpleConfSpace confSpace, Sequence seq, Parallelism parallelism, EnergyPartition epart, double epsilon, Result result) {

		// make the energy calculator
		var ecalc = new EnergyCalculator.Builder(new ForcefieldParams())
			.setParallelism(parallelism)
			.build();
		var confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
			.setEnergyPartition(epart)
			.build();

		// compute an emat
		result.ematStopwatch = new Stopwatch().start();
		var emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
			.build()
			.calcEnergyMatrix();
		result.ematStopwatch.stop();

		var rcs = seq.makeRCs(confSpace);
		gradientDescent(confEcalc, rcs, emat, ecalc.tasks, epsilon, result);
	}

	private static void gradientDescentCompiled(ConfSpace confSpace, Sequence seq, Parallelism parallelism, PosInterDist posInterDist, boolean includeStaticStatic, double epsilon, Result result) {

		try (var tasks = new ThreadPoolTaskExecutor()) {
			tasks.start(parallelism.numThreads);

			// make the energy calculator
			var ecalc = new CPUConfEnergyCalculator(confSpace);

			// compute an emat
			result.ematStopwatch = new Stopwatch().start();
			var emat = new EmatCalculator.Builder(ecalc)
				.setIncludeStaticStatic(includeStaticStatic)
				.setPosInterDist(posInterDist)
				.build()
				.calc();
			result.ematStopwatch.stop();

			var adapter = new ConfEnergyCalculatorAdapter.Builder(ecalc, tasks)
				.setIncludeStaticStatic(includeStaticStatic)
				.setPosInterDist(posInterDist)
				.build();
			var rcs = seq.makeRCs(confSpace);
			gradientDescent(adapter, rcs, emat, tasks, epsilon, result);
		}
	}

	private static void gradientDescent(ConfEnergyCalculator confEcalc, RCs rcs, EnergyMatrix emat, TaskExecutor tasks, double epsilon, Result result) {

		try (var ctx = tasks.contextGroup()) {

			var pfunc = new GradientDescentPfunc(
				confEcalc,
				new ConfAStarTree.Builder(emat, rcs)
					.setTraditional()
					.build(),
				new ConfAStarTree.Builder(emat, rcs)
					.setTraditional()
					.build(),
				rcs.getNumConformations()
			);
			pfunc.init(epsilon);
			pfunc.setInstanceId(0);
			pfunc.putTaskContexts(ctx);
			pfunc.setReportProgress(true);
			pfunc.setPreciseBcalc(true);

			result.pfuncStopwatch = new Stopwatch().start();
			pfunc.compute();
			result.pfuncStopwatch.stop();

			// finally, calculate the free energy
			var values = pfunc.getValues();
			result.freeEnergy = new FreeEnergyCalculator().calc(new MathTools.BigDecimalBounds(
				values.calcLowerBound(),
				values.calcUpperBound()
			));
		}
	}

	private static Result coffee(ConfSpace confSpace, Sequence seq, Parallelism parallelism, Bounds bounds, boolean includeStaticStatic, double gWidthMax) {

		var msConfSpace = new MultiStateConfSpace.Builder("complex", confSpace)
			.build();
		var state = msConfSpace.getState("complex");

		Coffee coffee = new Coffee.Builder(msConfSpace)
			.setParallelism(parallelism)
			.setStaticStatic(includeStaticStatic)
			//.setConditions(BoltzmannCalculator.Conditions.Body)
			//.setConditions(BoltzmannCalculator.Conditions.Room)
			.configEachState(config -> {
				config.ecalc = new CPUConfEnergyCalculator(confSpace);
				config.posInterGen = new PosInterGen(bounds.posInterDist, null);
			})
			.build();

		var director = new PfuncDirector.Builder(msConfSpace, state, seq)
			.setGWidthMax(gWidthMax)
			.setTiming(PfuncDirector.Timing.Precise)
			.build();

		var result = new Result();
		result.ematStopwatch = new Stopwatch().start();
		coffee.run(new Coffee.Director() {

			@Override
			public int numBestConfs() {
				return director.numBestConfs();
			}

			@Override
			public void init(Directions directions, NodeProcessor processor) {
				result.ematStopwatch.stop();
				director.init(directions, processor);
			}

			@Override
			public void direct(Directions directions, NodeProcessor processor) {
				result.pfuncStopwatch = new Stopwatch().start();
				director.direct(directions, processor);
				result.pfuncStopwatch.stop();
			}
		});
		result.freeEnergy = director.getFreeEnergy();

		return result;
	}
}
