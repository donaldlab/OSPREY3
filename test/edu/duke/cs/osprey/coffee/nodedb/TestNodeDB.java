package edu.duke.cs.osprey.coffee.nodedb;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.coffee.ClusterMember;
import edu.duke.cs.osprey.coffee.TestCoffee;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.MathTools;
import org.junit.Test;
import org.slf4j.bridge.SLF4JBridgeHandler;

import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.LongStream;


public class TestNodeDB {

	static {

		// configure hazelcast logging
		SLF4JBridgeHandler.removeHandlersForRootLogger();
		SLF4JBridgeHandler.install();
	}

	private static final long MiB = 1024*1024;

	private static void withMemNodeDB(MultiStateConfSpace confSpace, long dbBytes, Consumer<NodeDB> block) {
		withMemNodeDBs(confSpace, dbBytes, 1, block);
	}

	private static void withMemNodeDBs(MultiStateConfSpace confSpace, long dbBytes, int numMembers, Consumer<NodeDB> block) {
		var exceptions = ClusterMember.launchPseudoCluster(numMembers, cluster -> {
			try (var member = new ClusterMember(cluster)) {

				// make the node database
				try (var nodedb = new NodeDB.Builder(confSpace, member)
					.setMem(dbBytes)
					.build()
				) {

					// wait for all the database instances to be ready
					member.barrier(1, TimeUnit.MINUTES);

					block.accept(nodedb);
				}
			}
		});
		if (!exceptions.isEmpty()) {
			fail("Cluster threads encountered exceptions");
		}
	}

	@Test
	public void add1Poll1Local() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		withMemNodeDB(confSpace, MiB, nodedb -> {

			for (var state : confSpace.states) {

				// add a node to the local store
				var node = new NodeIndex.Node(
					state.index,
					Conf.make(state.confSpace),
					new BigExp(1.24, 6),
					new BigExp(6.8, 3)
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
		withMemNodeDB(confSpace, MiB, nodedb -> {

			var state = confSpace.states.get(0);

			// count dropped nodes
			var numDropped = new AtomicLong(0);
			nodedb.dropHandler = nodes -> numDropped.addAndGet(nodes.count());

			// add a bunch of random nodes
			Random rand = new Random(12345);
			for (int i=0; i<600_000; i++) {
				nodedb.addLocal(new NodeIndex.Node(
					state.index,
					Conf.make(state.confSpace),
					new BigExp(rand.nextDouble(), rand.nextInt()),
					new BigExp(rand.nextDouble(), rand.nextInt())
				));
			}

			// most of the nodes should get dropped
			assertThat(nodedb.size(state.index), greaterThan(10_000L));
			assertThat(nodedb.size(state.index), lessThan(600_000L));
			assertThat(numDropped.get(), greaterThan(400_000L));
			assertThat(numDropped.get(), lessThan(600_000L));
		});
	}

	private void addLotsLocalRemoveAll(int numThreads) {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		withMemNodeDB(confSpace, 10*MiB, nodedb -> {

			var state = confSpace.states.get(0);

			// count dropped nodes
			var numDropped = new AtomicLong(0);
			nodedb.dropHandler = nodes -> numDropped.addAndGet(nodes.count());

			// generate a bunch of random nodes
			var rand = new Random(12345);
			List<NodeIndex.Node> allNodes = LongStream.range(0, 100_000L)
				.mapToObj(i -> new NodeIndex.Node(
					state.index,
					Conf.make(state.confSpace),
					new BigExp(rand.nextDouble(), rand.nextInt()),
					new BigExp(rand.nextDouble(), rand.nextInt())
				))
				.collect(Collectors.toList());

			// add the nodes from a bunch of threads
			List<Thread> threads = IntStream.range(0, numThreads)
				.mapToObj(t -> new Thread(() -> {

					// add a slice of the nodes to the db
					int size = MathTools.divUp(allNodes.size(), numThreads);
					int start = t*size;
					int stop = Math.min(start + size, allNodes.size());
					for (var node : allNodes.subList(start, stop)) {
						nodedb.addLocal(node);
					}
				}))
				.collect(Collectors.toList());
			threads.forEach(t -> t.start());
			threads.forEach(t -> {
				try {
					t.join();
				} catch (InterruptedException ex) {
					throw new RuntimeException(ex);
				}
			});

			assertThat(numDropped.get(), is(0L));

			// sort them
			Comparator<NodeIndex.Node> comparator = Comparator.comparing(node -> node.score);
			var sortedNodes = new TreeSet<>(comparator);
			sortedNodes.addAll(allNodes);

			// remove all the nodes, check the scores
			assertThat(nodedb.size(state.index), is((long)allNodes.size()));
			for (int i=0; i<allNodes.size(); i++) {
				assertThat("" + i, nodedb.removeHigh(state.index).score, is(sortedNodes.pollLast().score));
			}
		});
	}
	@Test public void addLotsLocalRemoveAll_1() { addLotsLocalRemoveAll(1); }
	@Test public void addLotsLocalRemoveAll_2() { addLotsLocalRemoveAll(2); }
	@Test public void addLotsLocalRemoveAll_4() { addLotsLocalRemoveAll(4); }

	@Test
	public void addLocalRemoveHigh() {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();

		// make a random node
		var state = confSpace.states.get(0);
		Random rand = new Random(12345);
		var node = new NodeIndex.Node(
			state.index,
			Conf.make(state.confSpace),
			new BigExp(rand.nextDouble(), rand.nextInt()),
			new BigExp(rand.nextDouble(), rand.nextInt())
		);

		withMemNodeDBs(confSpace, MiB, 2, nodedb -> {

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

		withMemNodeDBs(confSpace, MiB, 2, nodedb -> {

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
						new BigExp(rand.nextDouble(), rand.nextInt()),
						new BigExp(rand.nextDouble(), rand.nextInt())
					));
				}
			}

			// wait for the node add to finish
			nodedb.member.barrier(10, TimeUnit.SECONDS);
			if (nodedb.member.id() != 0) {
				nodedb.member.waitForOperationsQuiet(1, TimeUnit.SECONDS, 5, TimeUnit.SECONDS);
			}
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

		withMemNodeDBs(confSpace, MiB, 4, nodedb -> {

			var state = confSpace.states.get(0);

			// add enough random nodes on each member to fill all the space
			// each member can hold ~55k nodes
			Random rand = new Random(12345);
			for (int i=0; i<60_000; i++) {
				nodedb.add(new NodeIndex.Node(
					state.index,
					Conf.make(state.confSpace),
					new BigExp(rand.nextDouble(), rand.nextInt()),
					new BigExp(rand.nextDouble(), rand.nextInt())
				));
			}

			// wait for the node add to finish
			nodedb.member.barrier(10, TimeUnit.SECONDS);
			nodedb.member.waitForOperationsQuiet(1, TimeUnit.SECONDS, 5, TimeUnit.SECONDS);

			// no one should have much free space left
			assertThat(nodedb.freeSpaceLocal(state.index), lessThan(nodedb.nodesPerBlock(state.index)));

			// wait for the query to finish
			nodedb.member.barrier(2, TimeUnit.SECONDS);
		});
	}
}
