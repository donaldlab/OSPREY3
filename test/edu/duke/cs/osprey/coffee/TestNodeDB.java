package edu.duke.cs.osprey.coffee;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.tools.BigExp;
import org.junit.Test;
import org.slf4j.bridge.SLF4JBridgeHandler;

import java.util.Random;
import java.util.function.Consumer;


public class TestNodeDB {

	static {

		// configure hazelcast logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();
	}

	private static final long MiB = 1024*1024;

	private static void withLocalMemNodeDB(MultiStateConfSpace confSpace, long dbBytes, Consumer<NodeDB> block) {
		try (var member = new ClusterMember(new Cluster("NodeDB", "job", 0, 1))) {
			NodeDB nodedb = new NodeDB(confSpace, member, null, 0, dbBytes);
			block.accept(nodedb);
		}
	}

	@Test
	public void add1Poll1Local() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		withLocalMemNodeDB(confSpace, 7*MiB, nodedb -> {

			for (var state : confSpace.states) {

				// add a node to the local store
				var node = new NodeIndex.Node(
					state.index,
					Conf.make(state.confSpace),
					new BigExp(1.24, 6)
				);
				nodedb.addLocal(node);
				assertThat(nodedb.size(state.index), is(1L));

				// remove the node
				var node2 = nodedb.pollHighestLocal(state.index);
				assertThat(node2, is(node));
				assertThat(nodedb.size(state.index), is(0L));
			}
		});
	}

	@Test
	public void fillStateLocal() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		withLocalMemNodeDB(confSpace, 7*MiB, nodedb -> {

			var state = confSpace.states.get(0);

			// add a bunch of random nodes
			Random rand = new Random(12345);
			for (int i=0; i<600_000; i++) {
				nodedb.addLocal(new NodeIndex.Node(
					state.index,
					Conf.make(state.confSpace),
					new BigExp(rand.nextDouble(), rand.nextInt())
				));
			}

			assertThat(nodedb.size(state.index), lessThan(600_000L));
		});
	}
}
