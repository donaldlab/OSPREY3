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
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
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
		//Total Z upper bound reduction through minimizations: 4.322064e+96
		//GD     classic   noSS   emat   17396 ms ( 17.40 s)   pfunc   77321 ms (  1.29 m)   G [-129.9658,-129.2924]  w =  0.6734

		//benchmark("GD     compiled  noSS", () -> gradientDescent(complexCompiled, seqCompiled, oneCpu, bounds, noStaticStatic, epsilon));
		//[14    40    28    37    28    13    13   ] scores:   18634, confs: 168, score: -101.009497, energy:  -99.118078, bounds:[   75.762014,   76.253323] (log10p1), delta:0.677380, time:   37.06 s, heapMem:3.5% of 1.9 GiB, extMem:0 B
		//GD     compiled  noSS   emat    9717 ms (  9.72 s)   pfunc   37045 ms ( 37.05 s)   G [-104.1276,-103.4567]  w =  0.6709
		// compiled is 2.09x faster than classic! =D

		//benchmark("GD     compiled yesSS", () -> gradientDescent(complexCompiled, seqCompiled, oneCpu, bounds, yesStaticStatic, epsilon));
		//[14    40    5     42    28    13    13   ] scores:   17994, confs: 167, score:-1557.838995, energy:-1553.935900, bounds:[ 1142.602893, 1143.096674] (log10p1), delta:0.679211, time:   36.57 s, heapMem:3.5% of 1.9 GiB, extMem:0 B
		//GD     compiled yesSS   emat    9490 ms (  9.49 s)   pfunc   36588 ms ( 36.59 s)   G [-1560.9539,-1560.2796]  w =  0.6743

		//benchmark("COFFEE compiled  noSS", () -> coffee(complexCompiled, seqCompiled, oneCpu, bounds, noStaticStatic, gWidthMax));
		//COFFEE-0: 	G [-104.119,-103.456]   width 0.663883 of 0.000000   confs       224   avgap 2.78   nodedb  22.7%   rr Infinity   time 50.08 s
		//COFFEE compiled  noSS   emat   17954 ms ( 17.95 s)   pfunc   50080 ms ( 50.08 s)   G [-104.1195,-103.4556]  w =  0.6639
		// COFFEE is 0.74x faster than gradient-descent! ;_;

		//benchmark("COFFEE compiled yesSS", () -> coffee(complexCompiled, seqCompiled, oneCpu, bounds, yesStaticStatic, gWidthMax));
		//COFFEE-0: 	G [-1560.929,-1560.271]   width 0.658356 of 0.000000   confs       258   avgap 2.82   nodedb  25.8%   rr Infinity   time 58.46 s
		//COFFEE compiled yesSS   emat   18845 ms ( 18.85 s)   pfunc   58459 ms ( 58.46 s)   G [-1560.9294,-1560.2710]  w =  0.6584

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
