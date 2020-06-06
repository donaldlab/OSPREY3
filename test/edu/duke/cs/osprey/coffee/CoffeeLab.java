package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.coffee.drivers.AffinityDriver;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.confspace.compiled.TestConfSpace;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import org.slf4j.bridge.SLF4JBridgeHandler;


public class CoffeeLab {

	public static void main(String[] args) {

		// configure hazelcast logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();

		ClusterMember.launchPseudoCluster(2, cluster -> run(cluster, Parallelism.makeCpu(2)));
	}

	private static void run(Cluster cluster, Parallelism parallelism) {

		// load a multi-state conf space
		var confSpaces = TestConfSpace.Design2RL0Interface7Mut.makeCompiled();
		MultiStateConfSpace confSpace = new MultiStateConfSpace
			.Builder("complex", confSpaces.complex)
			.addMutableState("design", confSpaces.chainA)
			.addMutableState("target", confSpaces.chainB)
			.build();

		var posInterDist = PosInterDist.DesmetEtAl1992;
		//var posInterDist = PosInterDist.TighterBounds;

		Coffee coffee = new Coffee.Builder(confSpace)
			.setCluster(cluster)
			.setParallelism(parallelism)
			.configEachState(config -> {
				config.ecalc = new CPUConfEnergyCalculator((ConfSpace)config.state.confSpace);
				config.posInterGen = new PosInterGen(posInterDist, null);
			})
			.build();

		var driver = new AffinityDriver(confSpace, "complex", "design", "target", 5);
		coffee.run(driver);
	}
}
