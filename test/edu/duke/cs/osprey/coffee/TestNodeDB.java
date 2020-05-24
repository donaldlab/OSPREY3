package edu.duke.cs.osprey.coffee;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.tools.BigExp;
import org.junit.Test;
import org.slf4j.bridge.SLF4JBridgeHandler;


public class TestNodeDB {

	static {

		// configure hazelcast logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();
	}

	@Test
	public void test() {

		// load a multi-state conf space
		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();

		final int statei = 0;
		try (var member = new ClusterMember(new Cluster("NodeDB", "job", 0, 1))) {
			NodeDB nodedb = new NodeDB(confSpace, member, null);

			nodedb.init(statei, Conf.make(confSpace.states.get(statei).confSpace), new BigExp(1.24, 6));
			// TODO
		}
	}
}
