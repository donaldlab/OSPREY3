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

		// single-threaded benchmarks

		benchmark("GD     classic   noSS", () -> gradientDescent(complexClassic, seqClassic, oneCpu, bounds, noStaticStatic, epsilon));
		//[28    4     0     28    13    28    10   ] scores:   20579, confs: 262, score: -126.800794, energy: -120.704369, bounds:[   94.681678,   95.174802] (log10p1), delta:0.678726, time:    1.52 m, heapMem:4.6% of 1.9 GiB, extMem:0 B
		//GD     classic   noSS   emat   20008 ms ( 20.01 s)   pfunc   91476 ms (  1.52 m)   G [-129.9658,-129.2924]  w =  0.6734

		benchmark("GD     compiled  noSS", () -> gradientDescent(complexCompiled, seqCompiled, oneCpu, bounds, noStaticStatic, epsilon));
		//[14    40    28    37    28    13    13   ] scores:   19469, confs: 168, score: -101.009497, energy:  -99.118078, bounds:[   75.762014,   76.253304] (log10p1), delta:0.677366, time:   48.90 s, heapMem:4.8% of 1.9 GiB, extMem:0 B
		//GD     compiled  noSS   emat   11748 ms ( 11.75 s)   pfunc   48934 ms ( 48.93 s)   G [-104.1276,-103.4567]  w =  0.6709
		// compiled is 1.86x faster than classic! =D

		benchmark("GD     compiled yesSS", () -> gradientDescent(complexCompiled, seqCompiled, oneCpu, bounds, yesStaticStatic, epsilon));
		//[14    40    28    37    28    13    13   ] scores:   19821, confs: 168, score:-1557.832821, energy:-1555.941402, bounds:[ 1142.603182, 1143.094467] (log10p1), delta:0.677362, time:   47.09 s, heapMem:5.0% of 1.9 GiB, extMem:0 B
		//GD     compiled yesSS   emat   13296 ms ( 13.30 s)   pfunc   47085 ms ( 47.09 s)   G [-1560.9509,-1560.2800]  w =  0.6709

		benchmark("COFFEE compiled  noSS", () -> coffee(complexCompiled, seqCompiled, oneCpu, bounds, noStaticStatic, gWidthMax));
		//COFFEE-0: 	G [-104.102,-103.452]   width 0.649504 of 0.000000   confs       257   avgap 2.80   nodedb  25.8%   rr Infinity   time 1.20 m
		//COFFEE compiled  noSS   emat   25818 ms ( 25.82 s)   pfunc   72077 ms (  1.20 m)   G [-104.1015,-103.4520]  w =  0.6495
		// COFFEE is 0.68x faster than gradient-descent! ;_;

		benchmark("COFFEE compiled yesSS", () -> coffee(complexCompiled, seqCompiled, oneCpu, bounds, yesStaticStatic, gWidthMax));
		//COFFEE-0: 	G [-1560.942,-1560.274]   width 0.668469 of 0.000000   confs       246   avgap 2.73   nodedb  26.6%   rr Infinity   time 1.33 m
		//COFFEE compiled yesSS   emat   18489 ms ( 18.49 s)   pfunc   79949 ms (  1.33 m)   G [-1560.9421,-1560.2736]  w =  0.6685

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
			.configEachState(config -> {
				config.ecalc = new CPUConfEnergyCalculator(confSpace);
				config.posInterGen = new PosInterGen(bounds.posInterDist, null);
			})
			.build();

		var director = new PfuncDirector(msConfSpace, state, msConfSpace.seqSpace.makeWildTypeSequence(), gWidthMax, PfuncDirector.Timing.Precise);

		var result = new Result();
		result.ematStopwatch = new Stopwatch().start();
		coffee.run(new Coffee.Director() {

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
