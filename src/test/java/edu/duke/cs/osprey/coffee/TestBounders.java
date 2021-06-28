package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.coffee.nodedb.NodeTree;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.MathTools.Optimizer;
import org.junit.Test;

import java.math.BigDecimal;
import java.util.Arrays;

import static edu.duke.cs.osprey.tools.Log.log;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;
import static edu.duke.cs.osprey.TestBase.isRelatively;

/**
 * Test the factor bounds.
 *
 * We should show that
 *  1. The factor bounds are always tighter than the old bounds,
 *  2. The factor bounds are always true bounds on the partition function
 *      * This is trivial for leaf nodes, and not necessary (is exactly the same gscore implementation)
 *      * For internal nodes I am satisfied by just confirming that we are true for "stem" nodes (1 from leaf)
 */
public class TestBounders {

    static {
        Cluster.fixHazelcastLogging();
    }

    @Test
    public void pos2_desmet_noSS_noTrip() {
        testNodeFactorBoundsLessThanOld(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "target",
                PosInterDist.DesmetEtAl1992,
                false,
                null,
                Bounds.AccurateUpper
        );
        testNodeFactorBoundsTruth(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "target",
                PosInterDist.DesmetEtAl1992,
                false,
                null,
                Bounds.AccurateUpper
        );
    }
    @Test
    public void pos2_desmet_noSS_noTrip_lower() {
        testNodeLowerBoundsGreaterThanOld(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "target",
                PosInterDist.DesmetEtAl1992,
                false,
                null,
                Bounds.Lower
        );
        testNodeLowerBoundsTruth(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "target",
                PosInterDist.DesmetEtAl1992,
                false,
                null,
                Bounds.Lower
        );
    }
    @Test
    public void pos2_desmet_SS_noTrip() {
        testNodeFactorBoundsLessThanOld(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "target",
                PosInterDist.DesmetEtAl1992,
                true,
                null,
                Bounds.AccurateUpper
        );
        testNodeFactorBoundsTruth(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "target",
                PosInterDist.DesmetEtAl1992,
                true,
                null,
                Bounds.AccurateUpper
        );
    }

    @Test
    public void pos2_desmet_SS_noTrip_lower() {
        testNodeLowerBoundsGreaterThanOld(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "target",
                PosInterDist.DesmetEtAl1992,
                true,
                null,
                Bounds.Lower
        );
        testNodeLowerBoundsTruth(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "target",
                PosInterDist.DesmetEtAl1992,
                true,
                null,
                Bounds.Lower
        );
    }

    // with only 2 positions, triples won't help, but it shouldn't crash
    @Test
    public void pos2_tighter_SS_trip() {
        testNodeFactorBoundsLessThanOld(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "target",
                PosInterDist.TighterBounds,
                true,
                1.0,
                Bounds.AccurateUpper
        );
        testNodeFactorBoundsTruth(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "target",
                PosInterDist.TighterBounds,
                true,
                1.0,
                Bounds.AccurateUpper
        );
    }
    @Test
    public void pos2_tighter_SS_trip_lower() {
        testNodeLowerBoundsGreaterThanOld(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "target",
                PosInterDist.TighterBounds,
                true,
                1.0,
                Bounds.Lower
        );
        testNodeLowerBoundsTruth(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "target",
                PosInterDist.TighterBounds,
                true,
                1.0,
                Bounds.Lower
        );
    }

    // with 3 positions, pairs should be accurate
    @Test
    public void pos3_desmet_SS_noTrip() {
        testNodeFactorBoundsLessThanOld(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.DesmetEtAl1992,
                true,
                null,
                Bounds.AccurateUpper
        );
    }
    // fails for [5, -1, 10], for both factor bounds and non-factor bounds.
    // This is probably a pairwise or single bound trespass
    @Test
    public void pos3_desmet_SS_noTrip_trespass() {
        testNodeFactorBoundsTruth(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.DesmetEtAl1992,
                true,
                null,
                Bounds.AccurateUpper
        );
    }
    @Test
    public void pos3_desmet_SS_noTrip_lower() {
        testNodeLowerBoundsGreaterThanOld(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.DesmetEtAl1992,
                true,
                null,
                Bounds.Lower
        );
        testNodeLowerBoundsTruth(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.DesmetEtAl1992,
                true,
                null,
                Bounds.Lower
        );
    }
    @Test
    public void pos3_tighter_SS_noTrip() {
        testNodeFactorBoundsLessThanOld(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.TighterBounds,
                true,
                null,
                Bounds.AccurateUpper
        );
        testNodeFactorBoundsTruth(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.TighterBounds,
                true,
                null,
                Bounds.AccurateUpper
        );
    }
    @Test
    public void pos3_tighter_SS_noTrip_lower() {
        testNodeLowerBoundsGreaterThanOld(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.TighterBounds,
                true,
                null,
                Bounds.Lower
        );
        testNodeLowerBoundsTruth(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.TighterBounds,
                true,
                null,
                Bounds.Lower
        );
    }

    // with 3 positions, triples should be accurate
    @Test
    public void pos3_desmet_SS_trip() {
        testNodeFactorBoundsLessThanOld(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.DesmetEtAl1992,
                true,
                1.0,
                Bounds.AccurateUpper
        );
    }
    // fails for [5, -1, 10], for both factor bounds and non-factor bounds.
    // This is probably a pairwise or single bound trespass
    @Test
    public void pos3_desmet_SS_trip_trespass(){
        // fails for [5, -1, 10]
        testNodeFactorBoundsTruth(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.DesmetEtAl1992,
                true,
                1.0,
                Bounds.AccurateUpper
        );
    }
    @Test
    public void pos3_desmet_SS_trip_lower() {
        testNodeLowerBoundsGreaterThanOld(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.DesmetEtAl1992,
                true,
                1.0,
                Bounds.Lower
        );
        testNodeLowerBoundsTruth(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.DesmetEtAl1992,
                true,
                1.0,
                Bounds.Lower
        );
    }

    // with 3 positions and "tighter" bonds
    @Test
    public void pos3_tighter_SS_trip() {
        testNodeFactorBoundsLessThanOld(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.TighterBounds,
                true,
                1.0,
                Bounds.AccurateUpper
        );
        testNodeFactorBoundsTruth(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.TighterBounds,
                true,
                1.0,
                Bounds.AccurateUpper
        );
    }
    @Test
    public void pos3_tighter_SS_trip_lower() {
        testNodeLowerBoundsGreaterThanOld(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.TighterBounds,
                true,
                1.0,
                Bounds.Lower
        );
        testNodeLowerBoundsTruth(
                TestCoffee.affinity_6ov7_1mut2flex(),
                "complex",
                PosInterDist.TighterBounds,
                true,
                1.0,
                Bounds.Lower
        );
    }

    private void testNodeFactorBoundsLessThanOld(MultiStateConfSpace confSpace, String stateName, PosInterDist posInterDist, boolean includeStaticStatic, Double tripleThreshold, Bounds bounds) {

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
        var zmatLower = coffee.calcZMat(state.index, true);


        final double epsilon = 1e-1;
        final int numNodes = 50;

        // get some confs
        var tree = new NodeTree(confSpace.seqSpace.makeWildTypeSequence().makeRCs(state.confSpace));
        var nodes = coffee.findFirstNNodes(state.index, zmat, tree, numNodes);

        var stateInfo = new StateInfo(stateConfig);
        stateInfo.zmat = zmat;
        stateInfo.initBounder();
        var confIndex = stateInfo.makeConfIndex();
        var factorStateInfo = new StateInfo(stateConfig);
        factorStateInfo.zmat = zmat;
        factorStateInfo.setFactorBounder(true);
        factorStateInfo.initBounder();
        var factorConfIndex = factorStateInfo.makeConfIndex();

        for (int i=0; i<nodes.size();i++){
            var node = nodes.get(i);
            Conf.index(node.conf, confIndex);
            Conf.index(node.conf, factorConfIndex);
            BigExp zold = stateInfo.zSumUpper(confIndex, tree);
            BigExp zfactor = factorStateInfo.zSumUpper(factorConfIndex, tree);
            // look at accuracy of bounds
            bounds.check(
                    String.format("%d/%d  %s: zold=%s ? zfactor=%s   (%f ? %f)",
                            i, nodes.size(), Conf.toString(node.conf),
                            zold, zfactor,
                            coffee.bcalc.freeEnergyPrecise(zold), coffee.bcalc.freeEnergyPrecise(zfactor)
                    ),
                    zold, zfactor, epsilon
            );
        }
    }

    private void testNodeLowerBoundsGreaterThanOld(MultiStateConfSpace confSpace, String stateName, PosInterDist posInterDist, boolean includeStaticStatic, Double tripleThreshold, Bounds bounds) {

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
        var zmatLower = coffee.calcZMat(state.index, true);


        final double epsilon = 1e-1;
        final int numNodes = 50;

        // get some confs
        var tree = new NodeTree(confSpace.seqSpace.makeWildTypeSequence().makeRCs(state.confSpace));
        var nodes = coffee.findFirstNNodes(state.index, zmat, tree, numNodes);

        var stateInfo = new StateInfo(stateConfig);
        stateInfo.zmat = zmat;
        stateInfo.zmatLower = zmatLower;
        stateInfo.initBounder();
        var confIndex = stateInfo.makeConfIndex();
        var factorStateInfo = new StateInfo(stateConfig);
        factorStateInfo.zmat = zmat;
        factorStateInfo.zmatLower = zmatLower;
        factorStateInfo.setFactorBounder(true);
        factorStateInfo.initBounder();
        var factorConfIndex = factorStateInfo.makeConfIndex();

        for (int i=0; i<nodes.size();i++){
            var node = nodes.get(i);
            Conf.index(node.conf, confIndex);
            Conf.index(node.conf, factorConfIndex);
            BigExp zold = stateInfo.zSumLower(confIndex, tree);
            BigExp zfactor = factorStateInfo.zSumLower(factorConfIndex, tree);
            // look at accuracy of bounds
            bounds.check(
                    String.format("%d/%d  %s: zold=%s ? zfactor=%s   (%f ? %f)",
                            i, nodes.size(), Conf.toString(node.conf),
                            zold, zfactor,
                            coffee.bcalc.freeEnergyPrecise(zold), coffee.bcalc.freeEnergyPrecise(zfactor)
                    ),
                    zold, zfactor, epsilon
            );
        }
    }

    private void testNodeFactorBoundsTruth(MultiStateConfSpace confSpace, String stateName, PosInterDist posInterDist, boolean includeStaticStatic, Double tripleThreshold, Bounds bounds) {
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


        final double epsilon = 1e-1;
        final int numNodes = 50;

        // get some confs
        var tree = new NodeTree(confSpace.seqSpace.makeWildTypeSequence().makeRCs(state.confSpace));
        var nodes = coffee.findStemNodesAndChildren(state.index, zmat, tree, numNodes);

        var stateInfo = new StateInfo(stateConfig);
        stateInfo.zmat = zmat;
        stateInfo.setFactorBounder(true);
        stateInfo.initBounder();
        var confIndex = stateInfo.makeConfIndex();

        try (var ecalc = new CPUConfEnergyCalculator(stateConfig.confSpace)) {
            for (int i = 0; i < nodes.size(); i++) {
                var nodePair = nodes.get(i);
                var node = nodePair.getKey();
                var children = nodePair.getValue();
                Conf.index(node.conf, confIndex);
                BigExp bound = stateInfo.zSumUpper(confIndex, tree);

                //int[] debugConf = {5, -1, 10};

                // minimize the children
                BigExp childSum = new BigExp(0.0);
                for (int j = 0; j < children.size(); j++){
                    var child = children.get(j);
                    var ecoords = ecalc.minimize(child.conf, stateConfig.makeInters(child.conf, includeStaticStatic));

                    /*
                    if (Arrays.equals(node.conf,debugConf)) {
                        childSum.add(new BigExp(coffee.bcalc.calcPrecise(ecoords.energy)));
                        System.out.println(String.format("child %d/%d  %s: %s (%f)",
                                j, children.size(), Conf.toString(child.conf),
                                coffee.bcalc.calcPrecise(ecoords.energy), ecoords.energy
                        ));

                        for (int posi1=0; posi1<stateConfig.confSpace.numPos(); posi1++) {
                            int confi1 = child.conf[posi1];
                            var boundContribution = zmat.single(posi1, confi1);
                            var inters = stateConfig.posInterGen.single(stateConfig.confSpace, posi1, confi1);
                            var z = new BigExp(coffee.bcalc.calcPrecise(ecalc.calcEnergy(ecoords.coords, inters)));
                            if (boundContribution.lessThan(z, epsilon)) {
                                log("WARNING: single bound %d:%d is wrong in conf %s: b=%s >=? z=%s   (%f <=? %f)",
                                        posi1, confi1, Conf.toString(child.conf),
                                        boundContribution, z,
                                        coffee.bcalc.freeEnergyPrecise(boundContribution), coffee.bcalc.freeEnergyPrecise(z)
                                );
                            }
                        }

                        // look at accuracy of pair bounds
                        for (int posi1=0; posi1<stateConfig.confSpace.numPos(); posi1++) {
                            int confi1 = child.conf[posi1];
                            for (int posi2=0; posi2<posi1; posi2++) {
                                int confi2 = child.conf[posi2];
                                var boundContribution = zmat.pair(posi1, confi1, posi2, confi2);
                                var inters = stateConfig.posInterGen.pair(stateConfig.confSpace, posi1, confi1, posi2, confi2);
                                var z = new BigExp(coffee.bcalc.calcPrecise(ecalc.calcEnergy(ecoords.coords, inters)));
                                if (boundContribution.lessThan(z, epsilon)) {
                                    log("WARNING: pair bound %d:%d %d:%d is wrong in conf %s: b=%s >=? z=%s   (%f <=? %f)",
                                            posi1, confi1, posi2, confi2, Conf.toString(child.conf),
                                            boundContribution, z,
                                            coffee.bcalc.freeEnergyPrecise(boundContribution), coffee.bcalc.freeEnergyPrecise(z)
                                    );
                                }
                            }
                        }
                    }

                     */
                }

                bounds.check(
                        String.format("%d/%d  %s: bound=%s ? childSum=%s   (%f ? %f)",
                                i, nodes.size(), Conf.toString(node.conf),
                                bound, childSum,
                                coffee.bcalc.freeEnergyPrecise(bound), coffee.bcalc.freeEnergyPrecise(childSum)
                        ),
                        bound, childSum, epsilon
                );
            }
        }
    }

    private void testNodeLowerBoundsTruth(MultiStateConfSpace confSpace, String stateName, PosInterDist posInterDist, boolean includeStaticStatic, Double tripleThreshold, Bounds bounds) {
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
        var zmatLower = coffee.calcZMat(state.index, true);


        final double epsilon = 1e-1;
        final int numNodes = 50;

        // get some confs
        var tree = new NodeTree(confSpace.seqSpace.makeWildTypeSequence().makeRCs(state.confSpace));
        var nodes = coffee.findStemNodesAndChildren(state.index, zmat, tree, numNodes);

        var stateInfo = new StateInfo(stateConfig);
        stateInfo.zmat = zmat;
        stateInfo.zmatLower = zmatLower;
        stateInfo.setFactorBounder(true);
        stateInfo.initBounder();
        var confIndex = stateInfo.makeConfIndex();

        try (var ecalc = new CPUConfEnergyCalculator(stateConfig.confSpace)) {
            for (int i = 0; i < nodes.size(); i++) {
                var nodePair = nodes.get(i);
                var node = nodePair.getKey();
                var children = nodePair.getValue();
                Conf.index(node.conf, confIndex);
                BigExp bound = stateInfo.zSumLower(confIndex, tree);

                // minimize the children
                BigExp childSum = new BigExp(0.0);
                for (int j = 0; j < children.size(); j++){
                    var child = children.get(j);
                    var ecoords = ecalc.minimize(child.conf, stateConfig.makeInters(child.conf, includeStaticStatic));
                    childSum.add(new BigExp(coffee.bcalc.calcPrecise(ecoords.energy)));
                    /*
                    System.out.println(String.format("child %d/%d  %s: %s (%f)",
                            j, children.size(), Conf.toString(child.conf),
                            coffee.bcalc.calcPrecise(ecoords.energy), ecoords.energy
                            ));
                     */
                }

                bounds.check(
                        String.format("%d/%d  %s: bound=%s ? childSum=%s   (%f ? %f)",
                                i, nodes.size(), Conf.toString(node.conf),
                                bound, childSum,
                                coffee.bcalc.freeEnergyPrecise(bound), coffee.bcalc.freeEnergyPrecise(childSum)
                        ),
                        bound, childSum, epsilon
                );
            }
        }
    }

    private enum Bounds {

        Perfect{
            @Override
            public void check(String msg, BigExp bound, BigExp z, double epsilon) {
                assertThat(msg, bound, isRelatively(z, epsilon));
            }
        },

        AccurateUpper {
            @Override
            public void check(String msg, BigExp bound, BigExp z, double epsilon) {
                assertThat(msg, bound.greaterThanOrEqual(z, epsilon), is(true));
            }
        },

        Lower {
            @Override
            public void check(String msg, BigExp bound, BigExp z, double epsilon) {
                assertThat(msg, bound.lessThanOrEqual(z), is(true));
            }
        };

        public abstract void check(String msg, BigExp bound, BigExp z, double epsilon);
    }
}
