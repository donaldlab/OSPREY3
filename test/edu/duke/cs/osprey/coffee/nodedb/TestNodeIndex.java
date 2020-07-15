package edu.duke.cs.osprey.coffee.nodedb;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.coffee.TestCoffee;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.BigExp;
import org.junit.Test;

import java.io.File;
import java.util.Comparator;
import java.util.Random;
import java.util.TreeSet;
import java.util.stream.IntStream;


public class TestNodeIndex {

	private void empty(File file) {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		var state = confSpace.getState("complex");

		var store = new BlockStore(file, 1024*1024);
		var index = new NodeIndex(store, state);

		assertThat(index.size(), is(0L));
		assertThat(index.removeHighest(), is(nullValue()));
	}

	@Test
	public void empty_mem() {
		empty(null);
	}

	@Test
	public void empty_file() {
		try (var file = new TestBase.TempFile("node.index")) {
			empty(file);
		}
	}

	private void addOneRemoveOne(File file) {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		var state = confSpace.getState("complex");

		var store = new BlockStore(file, 1024*1024);
		var index = new NodeIndex(store, state);

		var node = new NodeIndex.Node(
			state.index,
			Conf.make(state.confSpace),
			new BigExp(1.0, 42),
			new BigExp(4.2, 10)
		);
		boolean wasAdded = index.add(node);
		assertThat(wasAdded, is(true));

		assertThat(index.size(), is(1L));
		assertThat(index.highestScore(), is(node.score));
		assertThat(index.removeHighest(), is(node));

		assertThat(index.size(), is(0L));
		assertThat(index.highestScore(), is(nullValue()));
		assertThat(index.removeHighest(), is(nullValue()));
	}

	@Test
	public void addOneRemoveOne_mem() {
		addOneRemoveOne(null);
	}

	@Test
	public void addOneRemoveOne_file() {
		try (var file = new TestBase.TempFile("node.index")) {
			addOneRemoveOne(file);
		}
	}

	private void addFiveRemoveAll(File file) {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		var state = confSpace.getState("complex");

		var store = new BlockStore(file, 1024*1024);
		var index = new NodeIndex(store, state);

		var nodes = new NodeIndex.Node[] {
			new NodeIndex.Node(state.index, Conf.make(state.confSpace), new BigExp(1.0, 1), new BigExp(4.2, 5)),
			new NodeIndex.Node(state.index, Conf.make(state.confSpace), new BigExp(1.0, 2), new BigExp(4.2, 6)),
			new NodeIndex.Node(state.index, Conf.make(state.confSpace), new BigExp(1.0, 3), new BigExp(4.2, 7)),
			new NodeIndex.Node(state.index, Conf.make(state.confSpace), new BigExp(1.0, 4), new BigExp(4.2, 8)),
			new NodeIndex.Node(state.index, Conf.make(state.confSpace), new BigExp(1.0, 5), new BigExp(4.2, 9))
		};

		// add in ascending order
		index.add(nodes[0]);
		index.add(nodes[1]);
		index.add(nodes[2]);
		index.add(nodes[3]);
		index.add(nodes[4]);
		assertThat(index.size(), is(5L));

		assertThat(index.highestScore(), is(nodes[4].score));
		assertThat(index.removeHighest(), is(nodes[4]));
		assertThat(index.highestScore(), is(nodes[3].score));
		assertThat(index.removeHighest(), is(nodes[3]));
		assertThat(index.highestScore(), is(nodes[2].score));
		assertThat(index.removeHighest(), is(nodes[2]));
		assertThat(index.highestScore(), is(nodes[1].score));
		assertThat(index.removeHighest(), is(nodes[1]));
		assertThat(index.highestScore(), is(nodes[0].score));
		assertThat(index.removeHighest(), is(nodes[0]));
		assertThat(index.size(), is(0L));
		assertThat(index.removeHighest(), is(nullValue()));

		// add in descending order
		index.add(nodes[4]);
		index.add(nodes[3]);
		index.add(nodes[2]);
		index.add(nodes[1]);
		index.add(nodes[0]);
		assertThat(index.size(), is(5L));

		assertThat(index.highestScore(), is(nodes[4].score));
		assertThat(index.removeHighest(), is(nodes[4]));
		assertThat(index.highestScore(), is(nodes[3].score));
		assertThat(index.removeHighest(), is(nodes[3]));
		assertThat(index.highestScore(), is(nodes[2].score));
		assertThat(index.removeHighest(), is(nodes[2]));
		assertThat(index.highestScore(), is(nodes[1].score));
		assertThat(index.removeHighest(), is(nodes[1]));
		assertThat(index.highestScore(), is(nodes[0].score));
		assertThat(index.removeHighest(), is(nodes[0]));
		assertThat(index.size(), is(0L));
		assertThat(index.removeHighest(), is(nullValue()));
	}

	@Test
	public void addFiveRemoveAll_mem() {
		addFiveRemoveAll(null);
	}

	@Test
	public void addFiveRemoveAll_file() {
		try (var file = new TestBase.TempFile("node.index")) {
			addFiveRemoveAll(file);
		}
	}

	private void addBlocksRemoveAll(File file) {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		var state = confSpace.getState("complex");

		var store = new BlockStore(file, 1024*1024);
		var index = new NodeIndex(store, state);

		// make enough node to fill more than one block
		var nodesPerBlock = index.nodesPerBlock();
		var nodes = IntStream.range(0, nodesPerBlock*2)
			.mapToObj(i -> new NodeIndex.Node(state.index, Conf.make(state.confSpace), new BigExp(4.2, 7), new BigExp(1.0, i)))
			.toArray(NodeIndex.Node[]::new);

		// add in ascending order
		for (int i=0; i<nodes.length; i++) {
			index.add(nodes[i]);
		}
		assertThat(index.size(), is((long)nodes.length));

		for (int i=0; i<nodes.length; i++) {
			assertThat("node " + i, index.highestScore(), is(nodes[nodes.length - i - 1].score));
			assertThat("node " + i, index.removeHighest(), is(nodes[nodes.length - i - 1]));
		}
		assertThat(index.size(), is(0L));

		// add in descending order
		for (int i=0; i<nodes.length; i++) {
			index.add(nodes[nodes.length - i - 1]);
		}
		assertThat(index.size(), is((long)nodes.length));

		for (int i=0; i<nodes.length; i++) {
			assertThat("node " + i, index.highestScore(), is(nodes[nodes.length - i - 1].score));
			assertThat("node " + i, index.removeHighest(), is(nodes[nodes.length - i - 1]));
		}
		assertThat(index.size(), is(0L));
	}

	@Test
	public void addBlocksRemoveAll_mem() {
		addBlocksRemoveAll(null);
	}

	@Test
	public void addBlocksRemoveAll_file() {
		try (var file = new TestBase.TempFile("node.index")) {
			addBlocksRemoveAll(file);
		}
	}


	private void addLots(File file) {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		var state = confSpace.getState("complex");

		var store = new BlockStore(file, 1024*1024);
		var index = new NodeIndex(store, state);

		// add a bunch of random nodes
		Random rand = new Random(12345);
		for (int i=0; i<100_000; i++) {
			index.add(new NodeIndex.Node(
				state.index,
				Conf.make(state.confSpace),
				new BigExp(rand.nextDouble(), rand.nextInt()),
				new BigExp(rand.nextDouble(), rand.nextInt())
			));
			index.dropped().clear();
		}

		// there's only room for ~55k nodes in the index
		assertThat(index.size(), lessThan(60_000L));
	}

	@Test
	public void addLots_mem() {
		addLots(null);
	}

	@Test
	public void addLots_file() {
		try (var file = new TestBase.TempFile("node.index")) {
			addLots(file);
		}
	}

	public void addLotsRemoveAll(File file) {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		var state = confSpace.getState("complex");

		var store = new BlockStore(file, 4*1024*1024);
		var index = new NodeIndex(store, state);

		TreeSet<NodeIndex.Node> sortedNodes = new TreeSet<>(Comparator.comparing(node -> node.score));
		final int numNodes = 100_000;

		// add a bunch of random nodes
		Random rand = new Random(12345);
		for (int i=0; i<numNodes; i++) {

			NodeIndex.Node node = new NodeIndex.Node(
				state.index,
				Conf.make(state.confSpace),
				new BigExp(rand.nextDouble(), rand.nextInt()),
				new BigExp(rand.nextDouble(), rand.nextInt())
			);
			index.add(node);

			assertThat(index.dropped().size(), is(0));

			sortedNodes.add(node);
		}

		// poll some nodes, check the scores
		long size = index.size();
		for (int i=0; i<numNodes; i++) {
			assertThat("" + i, index.size(), is((long)sortedNodes.size()));
			assertThat("" + i, index.highestScore(), is(sortedNodes.last().score));
			assertThat("" + i, index.removeHighest(), is(sortedNodes.pollLast()));
			assertThat(index.dropped().size(), is(0));
		}

		assertThat(index.size(), is(size - numNodes));
	}

	@Test
	public void addLotsRemoveAll_mem() {
		addLotsRemoveAll(null);
	}

	@Test
	public void addLotsRemoveAll_file() {
		try (var file = new TestBase.TempFile("node.index")) {
			addLotsRemoveAll(file);
		}
	}

	public void fillRemoveSome(File file) {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		var state = confSpace.getState("complex");

		var store = new BlockStore(file, 1024*1024);
		var index = new NodeIndex(store, state);

		TreeSet<NodeIndex.Node> highestNodes = new TreeSet<>(Comparator.comparing(node -> node.score));
		final int numNodes = 10;

		// add a bunch of random nodes, fill up all the space
		Random rand = new Random(12345);
		for (int i=0; i<100_000; i++) {

			NodeIndex.Node node = new NodeIndex.Node(
				state.index,
				Conf.make(state.confSpace),
				new BigExp(rand.nextDouble(), rand.nextInt()),
				new BigExp(rand.nextDouble(), rand.nextInt())
			);
			index.add(node);
			index.dropped().clear();

			highestNodes.add(node);
			while (highestNodes.size() > numNodes) {
				highestNodes.pollFirst();
			}
		}

		// poll some nodes, check the scores
		long size = index.size();
		for (int i=0; i<numNodes; i++) {
			assertThat("" + i, index.highestScore(), is(highestNodes.last().score));
			assertThat("" + i, index.removeHighest(), is(highestNodes.pollLast()));
		}

		assertThat(index.size(), is(size - numNodes));
	}

	@Test
	public void fillRemoveSome_mem() {
		fillRemoveSome(null);
	}

	@Test
	public void fillRemoveSome_file() {
		try (var file = new TestBase.TempFile("node.index")) {
			fillRemoveSome(file);
		}
	}
}
