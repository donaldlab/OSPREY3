package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.coffee.directors.AffinityDirector;
import edu.duke.cs.osprey.coffee.directors.Timing;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import org.slf4j.bridge.SLF4JBridgeHandler;

import java.io.File;
import java.util.concurrent.TimeUnit;


public class CoffeeLab {

	public static void main(String[] args) {

		// configure hazelcast logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();

		int numMembers = 1;
		var parallelism = Parallelism.makeCpu(4);
		ClusterMember.launchPseudoCluster(numMembers, cluster -> run(cluster, parallelism));
	}

	private static void run(Cluster cluster, Parallelism parallelism) {

		// load a multi-state conf space
		var confSpace = TestCoffee.affinity_6ov7_2mut4flex();

		var posInterDist = PosInterDist.DesmetEtAl1992;
		//var posInterDist = PosInterDist.TighterBounds;

		Coffee coffee = new Coffee.Builder(confSpace)
			.setCluster(cluster)
			.setParallelism(parallelism)
			//.setSeqDBMathContext(new MathContext(4096, RoundingMode.HALF_UP))
			//.setTripleCorrectionThreshold(10.0)
			.configEachState((config, ecalc) -> {
				config.posInterGen = new PosInterGen(posInterDist, null);
			})
			.build();

		var director = new AffinityDirector.Builder(confSpace, "complex", "design", "target")
			//.setK(49)
			//.setK(28)
			.setK(10)
			.setTiming(Timing.Precise)
			.setEnsembleTracking(5, new File("ensembles"))
			.setEnsembleMinUpdate(30, TimeUnit.SECONDS)
			.build();
		coffee.run(director);
	}
}
