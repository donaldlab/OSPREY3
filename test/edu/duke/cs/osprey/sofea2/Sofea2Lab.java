package edu.duke.cs.osprey.sofea2;

import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.confspace.compiled.TestConfSpace;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.ForkCluster;
import edu.duke.cs.osprey.parallelism.Parallelism;
import org.slf4j.bridge.SLF4JBridgeHandler;


public class Sofea2Lab {

	public static void main(String[] args) {

		// configure hazelcast logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();

		// run in a pseudo-cluster
		ForkCluster.run(2, true, Sofea2Lab.class, Sofea2Lab::run);
	}

	private static void run(Cluster cluster) {

		Parallelism parallelism = Parallelism.makeCpu(2);

		// load a multi-state conf space
		var confSpaces = TestConfSpace.Design2RL0Interface7Mut.makeCompiled();
		MultiStateConfSpace confSpace = new MultiStateConfSpace
			.Builder("complex", confSpaces.complex)
			.addMutableState("design", confSpaces.chainA)
			.addMutableState("target", confSpaces.chainB)
			.build();

		var posInterDist = PosInterDist.DesmetEtAl1992;
		//var posInterDist = PosInterDist.TighterBounds;

		Sofea2 sofea = new Sofea2.Builder(confSpace)
			.setCluster(cluster)
			.setParallelism(parallelism)
			.configEachState(config -> {
				config.ecalc = new CPUConfEnergyCalculator((ConfSpace)config.state.confSpace);
				config.posInterGen = new PosInterGen(posInterDist, null);
			})
			.build();

		sofea.refine();
	}
}
