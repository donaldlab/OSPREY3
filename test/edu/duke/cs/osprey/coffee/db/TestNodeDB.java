package edu.duke.cs.osprey.coffee.db;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.coffee.ClusterMember;
import edu.duke.cs.osprey.coffee.TestCoffee;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.BigExp;
import org.junit.Test;
import org.slf4j.bridge.SLF4JBridgeHandler;

import java.util.Random;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;


public class TestNodeDB {

	static {

		// configure hazelcast logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();
	}

	private static final long MiB = 1024*1024;

	private static void withLocalMemNodeDB(MultiStateConfSpace confSpace, long dbBytes, Consumer<NodeDB> block) {
		withLocalMemNodeDBs(confSpace, dbBytes, 1, block);
	}

	private static void withLocalMemNodeDBs(MultiStateConfSpace confSpace, long dbBytes, int numMembers, Consumer<NodeDB> block) {
		var exceptions = ClusterMember.launchPseudoCluster(numMembers, member -> {

			// make the node database
			try (var nodedb = new NodeDB.Builder(confSpace, member)
				.setMem(dbBytes)
				.build()
			) {

				// wait for all the database instances to be ready
				member.barrier(1, TimeUnit.MINUTES);

				block.accept(nodedb);
			}
		});
		if (!exceptions.isEmpty()) {
			for (var t : exceptions) {
				t.printStackTrace();
			}
			fail("Cluster threads encountered exceptions");
		}
	}

	@Test
	public void add1Poll1Local() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		withLocalMemNodeDB(confSpace, MiB, nodedb -> {

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
				var node2 = nodedb.removeHighestLocal(state.index);
				assertThat(node2, is(node));
				assertThat(nodedb.size(state.index), is(0L));
			}
		});
	}

	@Test
	public void fillLocal() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		withLocalMemNodeDB(confSpace, MiB, nodedb -> {

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

	@Test
	public void addLocalRemoveHigh() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();

		// make a random node
		var state = confSpace.states.get(0);
		Random rand = new Random(12345);
		var node = new NodeIndex.Node(
			state.index,
			Conf.make(state.confSpace),
			new BigExp(rand.nextDouble(), rand.nextInt())
		);

		withLocalMemNodeDBs(confSpace, MiB, 2, nodedb -> {

			// add the node to member 0
			if (nodedb.member.id() == 0) {
				nodedb.addLocal(node);

				// force a broadcast now
				nodedb.broadcast();
			}

			// wait for the node add to finish
			nodedb.member.barrier(2, TimeUnit.SECONDS);

			// query it from member 1
			if (nodedb.member.id() == 1) {
				var node2 = nodedb.removeHigh(state.index);
				assertThat(node2, is(node));
			}

			// wait for the query to finish
			nodedb.member.barrier(2, TimeUnit.SECONDS);
		});
	}

	@Test
	public void fillCluster2Unidirectional() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();

		withLocalMemNodeDBs(confSpace, MiB, 2, nodedb -> {

			var state = confSpace.states.get(0);

			// add the nodes on member 0
			if (nodedb.member.id() == 0) {

				// add enough random nodes to fill all the space
				// each member can hold ~55k nodes
				Random rand = new Random(12345);
				for (int i=0; i<120_000; i++) {
					nodedb.add(new NodeIndex.Node(
						state.index,
						Conf.make(state.confSpace),
						new BigExp(rand.nextDouble(), rand.nextInt())
					));
				}
			}

			// wait for the node add to finish
			nodedb.member.barrier(10, TimeUnit.SECONDS);

			// no one should have much free space left
			assertThat(nodedb.member.name, nodedb.freeSpaceLocal(state.index), lessThan(nodedb.nodesPerBlock(state.index)));

			// wait for the query to finish
			nodedb.member.barrier(2, TimeUnit.SECONDS);
		});
	}

	@Test
	public void fillCluster4Omnidirectional() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();

		withLocalMemNodeDBs(confSpace, MiB, 4, nodedb -> {

			var state = confSpace.states.get(0);

			// add enough random nodes on each member to fill all the space
			// each member can hold ~55k nodes
			Random rand = new Random(12345);
			for (int i=0; i<60_000; i++) {
				nodedb.add(new NodeIndex.Node(
					state.index,
					Conf.make(state.confSpace),
					new BigExp(rand.nextDouble(), rand.nextInt())
				));
			}

			// wait for the node add to finish
			nodedb.member.barrier(10, TimeUnit.SECONDS);

			// no one should have much free space left
			assertThat(nodedb.freeSpaceLocal(state.index), lessThan(nodedb.nodesPerBlock(state.index)));

			// wait for the query to finish
			nodedb.member.barrier(2, TimeUnit.SECONDS);
		});
	}
}
