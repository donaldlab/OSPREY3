package edu.duke.cs.osprey.coffee;

import static edu.duke.cs.osprey.tools.Log.log;
import static org.hamcrest.Matchers.*;
import static org.hamcrest.MatcherAssert.*;
import static edu.duke.cs.osprey.TestBase.isRelatively;

import edu.duke.cs.osprey.coffee.nodedb.NodeTree;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.BigExp;
import org.junit.jupiter.api.Test;
import org.slf4j.bridge.SLF4JBridgeHandler;
import edu.duke.cs.osprey.tools.MathTools.Optimizer;


public class TestLeafNodeBounds {

	static {
		Cluster.fixHazelcastLogging();
	}

	// pairwise bounds on 2 positions should be accurate
	@Test
	public void pos2_desmet_noSS_noTrip() {
		checkLeafNodeBounds(
				TestCoffee.affinity_6ov7_1mut2flex(),
				"target",
				PosInterDist.DesmetEtAl1992,
				false,
				null,
				Bounds.Accurate
		);
	}
	@Test
	public void pos2_desmet_noSS_noTrip_lower() {
		checkLeafNodeBoundsLower(
				TestCoffee.affinity_6ov7_1mut2flex(),
				"target",
				PosInterDist.DesmetEtAl1992,
				false,
				null,
				Bounds.LowerBound
		);
	}
	@Test
	public void pos2_desmet_SS_noTrip() {
		checkLeafNodeBounds(
			TestCoffee.affinity_6ov7_1mut2flex(),
			"target",
			PosInterDist.DesmetEtAl1992,
			true,
			null,
			Bounds.Accurate
		);
	}
	@Test
	public void pos2_desmet_SS_noTrip_lower() {
		checkLeafNodeBoundsLower(
				TestCoffee.affinity_6ov7_1mut2flex(),
				"target",
				PosInterDist.DesmetEtAl1992,
				true,
				null,
				Bounds.LowerBound
		);
	}

	// with the "tighter" bounds on 2 positions, the bounds should be perfectly tight
	@Test
	public void pos2_tighter_noSS_noTrip() {
		checkLeafNodeBounds(
			TestCoffee.affinity_6ov7_1mut2flex(),
			"target",
			PosInterDist.TighterBounds,
			false,
			null,
			Bounds.Perfect
		);
	}
	@Test
	public void pos2_tighter_noSS_noTrip_lower() {
		checkLeafNodeBoundsLower(
				TestCoffee.affinity_6ov7_1mut2flex(),
				"target",
				PosInterDist.TighterBounds,
				false,
				null,
				Bounds.LowerBound
		);
	}
	@Test
	public void pos2_tighter_SS_noTrip() {
		checkLeafNodeBounds(
			TestCoffee.affinity_6ov7_1mut2flex(),
			"target",
			PosInterDist.TighterBounds,
			true,
			null,
			Bounds.Perfect
		);
	}
	@Test
	public void pos2_tighter_SS_noTrip_lower() {
		checkLeafNodeBoundsLower(
				TestCoffee.affinity_6ov7_1mut2flex(),
				"target",
				PosInterDist.TighterBounds,
				true,
				null,
				Bounds.LowerBound
		);
	}

	// with only 2 positions, triples won't help, but it shouldn't crash
	@Test
	public void pos2_tighter_SS_trip() {
		checkLeafNodeBounds(
			TestCoffee.affinity_6ov7_1mut2flex(),
			"target",
			PosInterDist.TighterBounds,
			true,
			1.0,
			Bounds.Perfect
		);
	}
	@Test
	public void pos2_tighter_SS_trip_lower() {
		checkLeafNodeBoundsLower(
				TestCoffee.affinity_6ov7_1mut2flex(),
				"target",
				PosInterDist.TighterBounds,
				true,
				1.0,
				Bounds.LowerBound
		);
	}

	// with 3 positions, pairs should be accurate
	@Test
	public void pos3_desmet_SS_noTrip() {
		checkLeafNodeBounds(
			TestCoffee.affinity_6ov7_1mut2flex(),
			"complex",
			PosInterDist.DesmetEtAl1992,
			true,
			null,
			Bounds.Accurate
		);
	}
	@Test
	public void pos3_desmet_SS_noTrip_lower() {
		checkLeafNodeBoundsLower(
				TestCoffee.affinity_6ov7_1mut2flex(),
				"complex",
				PosInterDist.DesmetEtAl1992,
				true,
				null,
				Bounds.LowerBound
		);
	}
	@Test
	public void pos3_tighter_SS_noTrip() {
		checkLeafNodeBounds(
			TestCoffee.affinity_6ov7_1mut2flex(),
			"complex",
			PosInterDist.TighterBounds,
			true,
			null,
			Bounds.Accurate
		);
	}
	@Test
	public void pos3_tighter_SS_noTrip_lower() {
		checkLeafNodeBoundsLower(
				TestCoffee.affinity_6ov7_1mut2flex(),
				"complex",
				PosInterDist.TighterBounds,
				true,
				null,
				Bounds.LowerBound
		);
	}

	// with 3 positions, triples should be accurate
	@Test
	public void pos3_desmet_SS_trip() {
		checkLeafNodeBounds(
			TestCoffee.affinity_6ov7_1mut2flex(),
			"complex",
			PosInterDist.DesmetEtAl1992,
			true,
			1.0,
			Bounds.Accurate
		);
	}
	@Test
	public void pos3_desmet_SS_trip_lower() {
		checkLeafNodeBoundsLower(
				TestCoffee.affinity_6ov7_1mut2flex(),
				"complex",
				PosInterDist.DesmetEtAl1992,
				true,
				1.0,
				Bounds.LowerBound
		);
	}

	// with 3 positions and "tighter" bonds, triples should be perfectly tight bounds
	@Test
	public void pos3_tighter_SS_trip() {
		checkLeafNodeBounds(
			TestCoffee.affinity_6ov7_1mut2flex(),
			"complex",
			PosInterDist.TighterBounds,
			true,
			1.0,
			Bounds.Perfect
		);
	}
	@Test
	public void pos3_tighter_SS_trip_lower() {
		checkLeafNodeBoundsLower(
				TestCoffee.affinity_6ov7_1mut2flex(),
				"complex",
				PosInterDist.TighterBounds,
				true,
				1.0,
				Bounds.LowerBound
		);
	}


	private void checkLeafNodeBounds(MultiStateConfSpace confSpace, String stateName, PosInterDist posInterDist, boolean includeStaticStatic, Double tripleThreshold, Bounds bounds) {

		var coffee = new Coffee.Builder(confSpace)
			.setNodeDBMem(16*1024*1024) // 16 MiB should be enough space for these tiny tests
			.setParallelism(Parallelism.makeCpu(Parallelism.getMaxNumCPUs()))
			.configEachState(config -> {
				config.posInterGen = new PosInterGen(posInterDist, null);
			})
			.setStaticStatic(includeStaticStatic)
			.setTripleCorrectionThreshold(tripleThreshold)
			.build();

		var state = confSpace.getState(stateName);
		var stateConfig = coffee.stateConfigs[state.index];

		// compute a zmat
		var zmat = coffee.calcZMat(state.index);

		final double epsilon = 1e-1; // TODO: can get tighter epsilon here?
		final int numNodes = 50;

		// get some confs
		var tree = new NodeTree(confSpace.seqSpace.makeWildTypeSequence().makeRCs(state.confSpace));
		var nodes = coffee.findHighZNodes(state.index, zmat, tree, numNodes);

		try (var ecalc = new CPUConfEnergyCalculator(stateConfig.confSpace)) {

			for (int i=0; i<nodes.size(); i++) {
				var node = nodes.get(i);

				// minimize the conformation
				var ecoords = ecalc.minimize(node.conf, stateConfig.makeInters(node.conf, includeStaticStatic));

				// look at accuracy of single bounds
				for (int posi1=0; posi1<stateConfig.confSpace.numPos(); posi1++) {
					int confi1 = node.conf[posi1];
					var bound = zmat.single(posi1, confi1);
					var inters = stateConfig.posInterGen.single(stateConfig.confSpace, posi1, confi1);
					var z = new BigExp(coffee.bcalc.calcPrecise(ecalc.calcEnergy(ecoords.coords, inters)));
					if (bound.lessThan(z, epsilon)) {
						log("WARNING: single bound %d:%d is wrong in conf %s: b=%s >=? z=%s   (%f <=? %f)",
							posi1, confi1, Conf.toString(node.conf),
							bound, z,
							coffee.bcalc.freeEnergyPrecise(bound), coffee.bcalc.freeEnergyPrecise(z)
						);
					}
				}

				// look at accuracy of pair bounds
				for (int posi1=0; posi1<stateConfig.confSpace.numPos(); posi1++) {
					int confi1 = node.conf[posi1];
					for (int posi2=0; posi2<posi1; posi2++) {
						int confi2 = node.conf[posi2];
						var bound = zmat.pair(posi1, confi1, posi2, confi2);
						var inters = stateConfig.posInterGen.pair(stateConfig.confSpace, posi1, confi1, posi2, confi2);
						var z = new BigExp(coffee.bcalc.calcPrecise(ecalc.calcEnergy(ecoords.coords, inters)));
						if (bound.lessThan(z, epsilon)) {
							log("WARNING: pair bound %d:%d %d:%d is wrong in conf %s: b=%s >=? z=%s   (%f <=? %f)",
								posi1, confi1, posi2, confi2, Conf.toString(node.conf),
								bound, z,
								coffee.bcalc.freeEnergyPrecise(bound), coffee.bcalc.freeEnergyPrecise(z)
							);
						}
					}
				}

				// check the overall bound
				var z = new BigExp(coffee.bcalc.calcPrecise(ecoords.energy));
				bounds.check(
					String.format("%d/%d  %s: bound=%s ? z=%s   (%f ? %f)",
						i, nodes.size(), Conf.toString(node.conf),
						node.zSumUpper, z,
						coffee.bcalc.freeEnergyPrecise(node.zSumUpper), coffee.bcalc.freeEnergyPrecise(z)
					),
					node.zSumUpper, z, epsilon
				);
			}
		}
	}
	private void checkLeafNodeBoundsLower(MultiStateConfSpace confSpace, String stateName, PosInterDist posInterDist, boolean includeStaticStatic, Double tripleThreshold, Bounds bounds) {

		var coffee = new Coffee.Builder(confSpace)
				.setNodeDBMem(16*1024*1024) // 16 MiB should be enough space for these tiny tests
				.setParallelism(Parallelism.makeCpu(Parallelism.getMaxNumCPUs()))
				.configEachState(config -> {
					config.posInterGen = new PosInterGen(posInterDist, null);
				})
				.setStaticStatic(includeStaticStatic)
				.setTripleCorrectionThreshold(tripleThreshold)
				.build();

		var state = confSpace.getState(stateName);
		var stateConfig = coffee.stateConfigs[state.index];
		var stateInfo = new StateInfo(stateConfig);

		// compute a zmat
		var zmat = coffee.calcZMat(state.index, false);
		var zmatLower = coffee.calcZMat(state.index, true);
		stateInfo.zmat = zmat;
		stateInfo.zmatLower = zmatLower;

		// make the lowerBounder
		stateInfo.setFactorBounder(false);
		stateInfo.initBounder();

		// make the confindex for scoring
		var confIndex = stateInfo.makeConfIndex();

		final double epsilon = 1e-1; // TODO: can get tighter epsilon here?
		final int numNodes = 50;

		// get some confs
		var tree = new NodeTree(confSpace.seqSpace.makeWildTypeSequence().makeRCs(state.confSpace));
		var nodes = coffee.findHighZNodes(state.index, zmat, tree, numNodes);

		try (var ecalc = new CPUConfEnergyCalculator(stateConfig.confSpace)) {

			for (int i=0; i<nodes.size(); i++) {
				var node = nodes.get(i);

				// minimize the conformation
				var ecoords = ecalc.minimize(node.conf, stateConfig.makeInters(node.conf, includeStaticStatic));

				// look at accuracy of single bounds
				for (int posi1=0; posi1<stateConfig.confSpace.numPos(); posi1++) {
					int confi1 = node.conf[posi1];
					var bound = zmatLower.single(posi1, confi1);
					var inters = stateConfig.posInterGen.single(stateConfig.confSpace, posi1, confi1);
					var z = new BigExp(coffee.bcalc.calcPrecise(ecalc.calcEnergy(ecoords.coords, inters)));
					if (bound.greaterThan(z, epsilon)) {
						log("WARNING: single bound %d:%d is wrong in conf %s: b=%s <=? z=%s   (%f >=? %f)",
								posi1, confi1, Conf.toString(node.conf),
								bound, z,
								coffee.bcalc.freeEnergyPrecise(bound), coffee.bcalc.freeEnergyPrecise(z)
						);
					}
				}

				// look at accuracy of pair bounds
				for (int posi1=0; posi1<stateConfig.confSpace.numPos(); posi1++) {
					int confi1 = node.conf[posi1];
					for (int posi2=0; posi2<posi1; posi2++) {
						int confi2 = node.conf[posi2];
						var bound = zmatLower.pair(posi1, confi1, posi2, confi2);
						var inters = stateConfig.posInterGen.pair(stateConfig.confSpace, posi1, confi1, posi2, confi2);
						var z = new BigExp(coffee.bcalc.calcPrecise(ecalc.calcEnergy(ecoords.coords, inters)));
						if (bound.greaterThan(z, epsilon)) {
							log("WARNING: pair bound %d:%d %d:%d is wrong in conf %s: b=%s >=? z=%s   (%f <=? %f)",
									posi1, confi1, posi2, confi2, Conf.toString(node.conf),
									bound, z,
									coffee.bcalc.freeEnergyPrecise(bound), coffee.bcalc.freeEnergyPrecise(z)
							);
						}
					}
				}

				// check the overall bound
                Conf.index(node.conf, confIndex);
				var lowerBound = stateInfo.zSumLower(confIndex, tree);
				var z = new BigExp(coffee.bcalc.calcPrecise(ecoords.energy));
				bounds.check(
						String.format("%d/%d  %s: bound=%s ? z=%s   (%f ? %f)",
								i, nodes.size(), Conf.toString(node.conf),
								lowerBound, z,
								coffee.bcalc.freeEnergyPrecise(lowerBound), coffee.bcalc.freeEnergyPrecise(z)
						),
						lowerBound, z, epsilon
				);
			}
		}
	}

	private enum Bounds {

		Perfect {
			@Override
			public void check(String msg, BigExp bound, BigExp z, double epsilon) {
				assertThat(msg, bound, isRelatively(z, epsilon));
			}
		},

		Accurate {
			@Override
			public void check(String msg, BigExp bound, BigExp z, double epsilon) {
				assertThat(msg, bound.greaterThanOrEqual(z, epsilon), is(true));
			}
		},

		LowerBound {
			@Override
			public void check(String msg, BigExp bound, BigExp z, double epsilon) {
				assertThat(msg, bound.lessThanOrEqual(z), is(true));
			}
		};

		public abstract void check(String msg, BigExp bound, BigExp z, double epsilon);
	}
}
