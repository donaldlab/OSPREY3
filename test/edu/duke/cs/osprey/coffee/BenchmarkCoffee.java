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
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculatorAdapter;
import edu.duke.cs.osprey.energy.compiled.NativeConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gpu.Structs;
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

		//vsClassic();
		//affinity_6ov7_1mut6flex();
		affinity_6ov7_1mut11flex();
	}

	private static void affinity_6ov7_1mut11flex() {

		// load a complex state
		var complex = (ConfSpace)TestCoffee.affinity_6ov7_1mut11flex().getState("complex").confSpace;

		// set some settings
		var seq = complex.seqSpace().makeWildTypeSequence();
		var staticStatic = true;
		//var bounds = Bounds.Classic;
		var bounds = Bounds.Tighter;
		//var precision = Structs.Precision.Float32;
		var precision = Structs.Precision.Float64;

		double gWidthMax = 1.0;

		// benchmarks, on jerry4 (up to 48 threads, 24 cores, 4 Titan V GPUs)
		//var parallelism = Parallelism.makeCpu(48);
		//var parallelism = Parallelism.make(48, 1);
		var parallelism = Parallelism.make(48, 4);

		benchmark("COFFEE", () -> coffee(complex, seq, parallelism, bounds, staticStatic, gWidthMax, 1024, precision));

		// cpus = 48
		//COFFEE-0: 	G [-1382.335,-1381.336]   width 0.999081 of 0.000000   confs     35534   avgap 4.86   nodedb  10.4%   rr Infinity   time 4.69 m
		//              COFFEE   emat   13060 ms ( 13.06 s)   pfunc  281446 ms (  4.69 m)   G [-1382.3351,-1381.3361]  w =  0.9991
		// 126.3 confs/s

		// cpus = 48, precision = 32
		//COFFEE-0: 	G [-1382.253,-1381.257]   width 0.996068 of 0.000000   confs     29927   avgap 4.56   nodedb   7.7%   rr Infinity   time 2.71 m
		//              COFFEE   emat   10534 ms ( 10.53 s)   pfunc  162315 ms (  2.71 m)   G [-1382.2527,-1381.2566]  w =  0.9961
		// 184.4 confs/s

		// cpus = 48, gpus = 1
		//COFFEE-0: 	G [-1382.328,-1381.329]   width 0.998872 of 0.000000   confs     43769   avgap 4.91   nodedb   5.9%   rr Infinity   time 1.92 m
		//              COFFEE   emat   14172 ms ( 14.17 s)   pfunc  115357 ms (  1.92 m)   G [-1382.3282,-1381.3294]  w =  0.9989
		// 379.4 conf/s

		// cpus = 48, gpus = 4
		//COFFEE-0: 	G [-1382.318,-1381.330]   width 0.988346 of 0.000000   confs     59701   avgap 4.97   nodedb   4.0%   rr Infinity   time 52.29 s
		//              COFFEE   emat    9916 ms (  9.92 s)   pfunc   52294 ms ( 52.29 s)   G [-1382.3184,-1381.3301]  w =  0.9883
		// 1141.6 confs/s

		// cpus = 48, gpus = 4, precision = 32
		//COFFEE-0: 	G [-1382.274,-1381.278]   width 0.995364 of 0.000000   confs     48500   avgap 4.56   nodedb   3.3%   rr Infinity   time 48.25 s
		//              COFFEE   emat   12634 ms ( 12.63 s)   pfunc   48257 ms ( 48.26 s)   G [-1382.2738,-1381.2784]  w =  0.9954
		// 1005.0 confs/s

		// float32 isn't much faster... so minimizations aren't the bottleneck? or the node scoring heuristics need some improvement?
	}

	private static void affinity_6ov7_1mut6flex() {

		// load a complex state
		var complex = (ConfSpace)TestCoffee.affinity_6ov7_1mut6flex().getState("complex").confSpace;

		// set some settings
		var seq = complex.seqSpace().makeWildTypeSequence();
		var staticStatic = true;
		var bounds = Bounds.Classic;
		var precision = Structs.Precision.Float64;

		double epsilon = 0.68;
		double gWidthMax = 0.67; // amazingly corresponds to epsilon ~0.68

		// benchmarks, on jerry4 (up to 48 threads, 24 cores, 4 Titan V GPUs)
		//var parallelism = Parallelism.makeCpu(4);
		//var parallelism = Parallelism.makeCpu(48);
		//var parallelism = Parallelism.make(3, 1, 1);
		var parallelism = Parallelism.make(3*4, 4, 1);

		benchmark("GrdDsc", () -> gradientDescent(complex, seq, parallelism, bounds, staticStatic, epsilon));
		benchmark("COFFEE", () -> coffee(complex, seq, parallelism, bounds, staticStatic, gWidthMax, 2, precision));

		// cpus = 4
		//8     8     0     17    16    4    ] scores:  101383, confs:3418, score:-1376.578351, energy:-1373.434374, bounds:[ 1010.963826, 1011.457663] (log10p1), delta:0.679253, time:    1.75 m, heapMem:0.1% of 30.0 GiB, extMem:0 B
		//              GrdDsc   emat    2469 ms (  2.47 s)   pfunc  104917 ms (  1.75 m)   G [-1381.1945,-1380.5201]  w =  0.6744
		//COFFEE-0: 	G [-1381.182,-1380.520]   width 0.662268 of 0.000000   confs      3707   avgap 4.52   nodedb  94.5%   rr Infinity   time 1.96 m
		//              COFFEE   emat    5594 ms (  5.59 s)   pfunc  117770 ms (  1.96 m)   G [-1381.1824,-1380.5201]  w =  0.6623

		// cpus = 24
		//[7     8     0     4     11    0    ] scores:  131234, confs:3436, score:-1376.585883, energy:-1370.705288, bounds:[ 1010.963863, 1011.454391] (log10p1), delta:0.676800, time:   24.23 s, heapMem:0.1% of 30.0 GiB, extMem:0 B
		//              GrdDsc   emat    2252 ms (  2.25 s)   pfunc   24159 ms ( 24.16 s)   G [-1381.1900,-1380.5202]  w =  0.6698
		//COFFEE-0: 	G [-1381.181,-1380.520]   width 0.660512 of 0.000000   confs      3680   avgap 4.55   nodedb  93.0%   rr Infinity   time 25.05 s
		//              COFFEE   emat    4764 ms (  4.76 s)   pfunc   25055 ms ( 25.06 s)   G [-1381.1806,-1380.5201]  w =  0.6605

		// cpus = 48
		//[7     8     0     4     11    0    ] scores:  154288, confs:3459, score:-1376.585883, energy:-1370.705288, bounds:[ 1010.963883, 1011.450285] (log10p1), delta:0.673714, time:   19.99 s, heapMem:0.1% of 30.0 GiB, extMem:0 B
		//              GrdDsc   emat    2312 ms (  2.31 s)   pfunc   19865 ms ( 19.87 s)   G [-1381.1844,-1380.5202]  w =  0.6642
		//COFFEE-0: 	G [-1381.144,-1380.520]   width 0.623509 of 0.000000   confs      3775   avgap 4.61   nodedb  87.5%   rr Infinity   time 22.95 s
		//              COFFEE   emat    5330 ms (  5.33 s)   pfunc   22951 ms ( 22.95 s)   G [-1381.1439,-1380.5204]  w =  0.6235

		// gpus = 1
		//COFFEE-0: 	G [-1380.944,-1380.513]   width 0.431091 of 0.000000   confs      5376   avgap 4.78   nodedb  46.9%   rr Infinity   time 10.02 s
		//              COFFEE   emat    9889 ms (  9.89 s)   pfunc   10024 ms ( 10.02 s)   G [-1380.9436,-1380.5125]  w =  0.4311

		// gpus = 4
		//COFFEE-0: 	G [-1381.090,-1380.511]   width 0.578638 of 0.000000   confs      4106   avgap 4.71   nodedb  51.6%   rr Infinity   time 7.60 s
		//              COFFEE   emat   11785 ms ( 11.79 s)   pfunc    7607 ms (  7.61 s)   G [-1381.0901,-1380.5114]  w =  0.5786
	}

	private static void vsClassic() {

		// load a complex state
		var complexClassic = TestConfSpace.Design2RL0Interface7Mut.makeClassic().complex;
		var complexCompiled = TestConfSpace.Design2RL0Interface7Mut.makeCompiled().complex;

		var seqClassic = complexClassic.seqSpace.makeWildTypeSequence();
		var seqCompiled = complexCompiled.seqSpace.makeWildTypeSequence();

		var bounds = Bounds.Classic;
		var precision = Structs.Precision.Float64;

		var noStaticStatic = false;
		var yesStaticStatic = true;

		double epsilon = 0.68; // amazingly, also corresponds to roughly a 0.67 width in free energy
		//double epsilon = 0.018; // roughly equivalent to gWidthMax = 0.01
		double gWidthMax = 0.67;
		//double gWidthMax = 0.01;

		// benchmarks, on jerry4 (up to 48 threads, 24 cores, 4 Titan V GPUs)
		//var parallelism = Parallelism.makeCpu(1);
		//var parallelism = Parallelism.makeCpu(48);
		//var parallelism = Parallelism.make(48, 1, 64);
		var parallelism = Parallelism.make(48, 4, 64);

		//benchmark("GD     classic   noSS", () -> gradientDescent(complexClassic, seqClassic, parallelism, bounds, noStaticStatic, epsilon));

		// 1 thread
		//[28    4     0     28    13    28    10   ] scores:   16421, confs: 262, score: -126.800794, energy: -120.704369, bounds:[   94.681678,   95.174941] (log10p1), delta:0.678828, time:    1.32 m, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     classic   noSS   emat   19095 ms ( 19.10 s)   pfunc   79191 ms (  1.32 m)   G [-129.9660,-129.2924]  w =  0.6736

		// 2 threads
		//[28    13    0     28    9     28    6    ] scores:   19661, confs: 264, score: -126.799888, energy: -123.291188, bounds:[   94.681851,   95.170677] (log10p1), delta:0.675531, time:   41.12 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     classic   noSS   emat   10805 ms ( 10.81 s)   pfunc   41114 ms ( 41.11 s)   G [-129.9602,-129.2927]  w =  0.6675

		// 6 threads
		//[28    13    13    2     10    28    12   ] scores:   20212, confs: 268, score: -126.773173, energy: -124.492100, bounds:[   94.682038,   95.162465] (log10p1), delta:0.669194, time:   15.15 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     classic   noSS   emat    4287 ms (  4.29 s)   pfunc   15149 ms ( 15.15 s)   G [-129.9490,-129.2929]  w =  0.6560

		// 12 threads
		//[28    4     13    2     11    28    8    ] scores:   23299, confs: 274, score: -126.731286, energy: -124.574063, bounds:[   94.682310,   95.150630] (log10p1), delta:0.659842, time:    8.25 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     classic   noSS   emat    2988 ms (  2.99 s)   pfunc    8241 ms (  8.24 s)   G [-129.9328,-129.2933]  w =  0.6395

		// 24 threads
		//[28    4     0     2     11    28    8    ] scores:   23183, confs: 287, score: -126.687537, energy: -124.058369, bounds:[   94.682685,   95.125548] (log10p1), delta:0.639308, time:    5.96 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     classic   noSS   emat    2416 ms (  2.42 s)   pfunc    5954 ms (  5.95 s)   G [-129.8985,-129.2938]  w =  0.6048

		// 48 threads
		//[28    13    13    2     10    1     43   ] scores:   24442, confs: 311, score: -126.592165, energy: -122.776608, bounds:[   94.683647,   95.082771] (log10p1), delta:0.601089, time:    5.57 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     classic   noSS   emat    2305 ms (  2.31 s)   pfunc    5565 ms (  5.57 s)   G [-129.8401,-129.2951]  w =  0.5450

		//benchmark("GD     compiled yesSS", () -> gradientDescent(complexCompiled, seqCompiled, parallelism, bounds, yesStaticStatic, epsilon));

		// 1 thread
		//[28    38    28    39    28    13    9    ] scores:   12673, confs: 169, score:-1557.822975, energy:-1555.779307, bounds:[ 1142.603402, 1143.093190] (log10p1), delta:0.676249, time:   28.58 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     compiled yesSS   emat    8521 ms (  8.52 s)   pfunc   28465 ms ( 28.47 s)   G [-1560.9492,-1560.2803]  w =  0.6688

		// 2 threads
		//[14    38    28    41    28    13    9    ] scores:   10563, confs: 171, score:-1557.808425, energy:-1555.769221, bounds:[ 1142.603648, 1143.091199] (log10p1), delta:0.674576, time:   15.19 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     compiled yesSS   emat    8383 ms (  8.38 s)   pfunc   15112 ms ( 15.11 s)   G [-1560.9464,-1560.2807]  w =  0.6658

		// 6 threads
		//[14    38    28    41    28    13    9    ] scores:   11118, confs: 175, score:-1557.808425, energy:-1555.769221, bounds:[ 1142.604324, 1143.081718] (log10p1), delta:0.666876, time:    5.64 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     compiled yesSS   emat    8513 ms (  8.51 s)   pfunc    5568 ms (  5.57 s)   G [-1560.9335,-1560.2816]  w =  0.6519

		// 12 threads
		//[14    38    28    41    28    13    9    ] scores:   10639, confs: 180, score:-1557.808425, energy:-1555.769221, bounds:[ 1142.604855, 1143.071699] (log10p1), delta:0.658685, time:    3.21 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     compiled yesSS   emat    8522 ms (  8.52 s)   pfunc    3129 ms (  3.13 s)   G [-1560.9198,-1560.2823]  w =  0.6375

		// 24 threads
		//[14    36    28    41    28    13    13   ] scores:   10021, confs: 193, score:-1557.750911, energy:-1555.812669, bounds:[ 1142.606972, 1143.046281] (log10p1), delta:0.636344, time:    2.28 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     compiled yesSS   emat    8438 ms (  8.44 s)   pfunc    2199 ms (  2.20 s)   G [-1560.8851,-1560.2852]  w =  0.5999

		// 48 threads
		//[14    38    28    39    28    13    9    ] scores:    8083, confs: 221, score:-1557.696283, energy:-1555.650422, bounds:[ 1142.610232, 1143.004953] (log10p1), delta:0.597024, time:    2.29 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     compiled yesSS   emat    8491 ms (  8.49 s)   pfunc    2217 ms (  2.22 s)   G [-1560.8287,-1560.2897]  w =  0.5390

		// with e <= 0.01

		// 48 threads
		//[28    38    28    42    15    13    7    ] scores:   55323, confs:1514, score:-1554.136302, energy:-1551.137430, bounds:[ 1142.622923, 1142.630201] (log10p1), delta:0.016618, time:   12.83 s, heapMem:0.2% of 30.0 GiB, extMem:0 B
		//GD     compiled yesSS   emat    8467 ms (  8.47 s)   pfunc   12756 ms ( 12.76 s)   G [-1560.3169,-1560.3070]  w =  0.0099


		benchmark("COFFEE compiled yesSS", () -> coffee(complexCompiled, seqCompiled, parallelism, bounds, yesStaticStatic, gWidthMax, 2, precision));

		// 1 thread
		//COFFEE-0: 	G [-1560.934,-1560.277]   width 0.656849 of 0.000000   confs       226   avgap 2.75   nodedb  21.1%   rr Infinity   time 39.02 s
		//COFFEE compiled yesSS   emat   20914 ms ( 20.91 s)   pfunc   39022 ms ( 39.02 s)   G [-1560.9342,-1560.2773]  w =  0.6568

		// 2 threads
		//COFFEE-0: 	G [-1560.939,-1560.278]   width 0.660169 of 0.000000   confs       231   avgap 2.75   nodedb  20.3%   rr Infinity   time 21.57 s
		//COFFEE compiled yesSS   emat   11044 ms ( 11.04 s)   pfunc   21576 ms ( 21.58 s)   G [-1560.9386,-1560.2784]  w =  0.6602

		// 6 threads
		//COFFEE-0: 	G [-1560.915,-1560.281]   width 0.634783 of 0.000000   confs       231   avgap 2.68   nodedb  18.8%   rr Infinity   time 8.25 s
		//COFFEE compiled yesSS   emat    7771 ms (  7.77 s)   pfunc    8251 ms (  8.25 s)   G [-1560.9153,-1560.2806]  w =  0.6348

		// 12 threads
		//COFFEE-0: 	G [-1560.886,-1560.282]   width 0.603408 of 0.000000   confs       245   avgap 2.72   nodedb  19.5%   rr Infinity   time 5.12 s
		//COFFEE compiled yesSS   emat    7133 ms (  7.13 s)   pfunc    5123 ms (  5.12 s)   G [-1560.8857,-1560.2823]  w =  0.6034

		// 24 threads
		//COFFEE-0: 	G [-1560.704,-1560.297]   width 0.406967 of 0.000000   confs       294   avgap 2.73   nodedb  20.3%   rr Infinity   time 3.78 s
		//COFFEE compiled yesSS   emat    7149 ms (  7.15 s)   pfunc    3787 ms (  3.79 s)   G [-1560.7038,-1560.2968]  w =  0.4070

		// 48 threads
		//COFFEE-0: 	G [-1560.920,-1560.287]   width 0.633140 of 0.000000   confs       210   avgap 2.64   nodedb  18.8%   rr Infinity   time 2.86 s
		//COFFEE compiled yesSS   emat    8689 ms (  8.69 s)   pfunc    2858 ms (  2.86 s)   G [-1560.9200,-1560.2868]  w =  0.6331

		// 1 GPU
		//COFFEE-0: 	G [-1560.566,-1560.298]   width 0.267274 of 0.000000   confs       908   avgap 3.00   nodedb  17.2%   rr Infinity   time 2.96 s
		//COFFEE compiled yesSS   emat    6443 ms (  6.44 s)   pfunc    2964 ms (  2.96 s)   G [-1560.5656,-1560.2983]  w =  0.2673

		// with w <= 0.01

		// 48 threads
		//COFFEE-0: 	G [-1560.316,-1560.307]   width 0.008698 of 0.000000   confs      1649   avgap 2.93   nodedb  44.5%   rr Infinity   time 15.80 s
		//COFFEE compiled yesSS   emat    8848 ms (  8.85 s)   pfunc   15804 ms ( 15.80 s)   G [-1560.3157,-1560.3070]  w =  0.0087

		// 1 GPU
		//COFFEE-0: 	G [-1560.301,-1560.300]   width 0.001530 of 0.000000   confs      3242   avgap 2.93   nodedb  38.3%   rr Infinity   time 9.61 s
		//COFFEE compiled yesSS   emat   11113 ms ( 11.11 s)   pfunc    9612 ms (  9.61 s)   G [-1560.3012,-1560.2996]  w =  0.0015

		// 4 GPUs
		//COFFEE-0: 	G [-1560.300,-1560.300]   width 0.000627 of 0.000000   confs      5181   avgap 2.94   nodedb  39.1%   rr Infinity   time 4.63 s
		//COFFEE compiled yesSS   emat   12829 ms ( 12.83 s)   pfunc    4632 ms (  4.63 s)   G [-1560.3003,-1560.2997]  w =  0.0006
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

		// if the benchmark got skipped, don't report anything
		if (result.ematStopwatch == null || result.pfuncStopwatch == null) {
			return;
		}

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

		if (parallelism.numGpus > 0) {
			log("GrdDsc isn't configured to handle GPUs well, skipping benchmark");
			return;
		}

		try (var tasks = new ThreadPoolTaskExecutor()) {
			tasks.start(parallelism.numThreads);

			// make the energy calculator
			var ecalc = new NativeConfEnergyCalculator(confSpace, Structs.Precision.Float64);

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

	private static Result coffee(ConfSpace confSpace, Sequence seq, Parallelism parallelism, Bounds bounds, boolean includeStaticStatic, double gWidthMax, int nodesMiB, Structs.Precision precision) {

		var msConfSpace = new MultiStateConfSpace.Builder("complex", confSpace)
			.build();
		var state = msConfSpace.getState("complex");

		Coffee coffee = new Coffee.Builder(msConfSpace)
			.setParallelism(parallelism)
			.setStaticStatic(includeStaticStatic)
			//.setConditions(BoltzmannCalculator.Conditions.Body)
			//.setConditions(BoltzmannCalculator.Conditions.Room)
			.setNodeDBMem(nodesMiB*1024*1024)
			.setPrecision(precision)
			.configEachState((config, ecalc) -> {
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
