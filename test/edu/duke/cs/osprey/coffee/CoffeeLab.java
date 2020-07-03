package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.coffee.directors.SequenceDirector;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import org.slf4j.bridge.SLF4JBridgeHandler;

import static edu.duke.cs.osprey.tools.Log.log;


public class CoffeeLab {

	public static void main(String[] args) {

		// configure hazelcast logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();

		ClusterMember.launchPseudoCluster(2, cluster -> run(cluster, Parallelism.makeCpu(2)));
	}

	private static void run(Cluster cluster, Parallelism parallelism) {

		// load a multi-state conf space
		var confSpace = TestCoffee.affinity_6ov7_1mut2flex();

		var posInterDist = PosInterDist.DesmetEtAl1992;
		//var posInterDist = PosInterDist.TighterBounds;

		Coffee coffee = new Coffee.Builder(confSpace)
			.setCluster(cluster)
			.setParallelism(parallelism)
			.configEachState((config, ecalc) -> {
				config.posInterGen = new PosInterGen(posInterDist, null);
			})
			.build();

		//var director = new AffinityDirector(confSpace, "complex", "design", "target", 5);
		var director = new SequenceDirector(confSpace, confSpace.seqSpace.makeWildTypeSequence(), 0.1);
		coffee.run(director);

		// show the results
		for (var state : confSpace.states) {
			var g = director.getFreeEnergy(state);
			log("%20s   G = %s  w=%.2f", state.name, g, g.size());
		}
	}
}
