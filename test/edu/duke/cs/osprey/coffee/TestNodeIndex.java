package edu.duke.cs.osprey.coffee;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.BigExp;
import org.junit.Test;

import java.io.File;
import java.util.Random;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Consumer;


public class TestNodeIndex {

	private void add(File file) {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		var state = confSpace.getState("complex");

		// track the number of evicted nodes
		AtomicLong numEvicted = new AtomicLong(0);
		Consumer<NodeIndex.Node> evictionListener = node -> numEvicted.incrementAndGet();

		var db = new FixedDB(file, 2*1024*1024);
		var index = new NodeIndex(db, "index", state, evictionListener);

		// add a bunch of random nodes
		Random rand = new Random(12345);
		for (int i=0; i<100_000; i++) {
			index.add(new NodeIndex.Node(
				state.index,
				Conf.make(state.confSpace),
				new BigExp(rand.nextDouble(), rand.nextInt())
			));
		}

		assertThat(numEvicted.get(), greaterThan(0L));

		// there's only room for ~30k nodes in the index
		assertThat(index.size(), lessThan(40_000L));
	}

	@Test
	public void add_mem() {
		add(null);
	}

	@Test
	public void add_file() {
		try (var file = new TestBase.TempFile("node.index")) {
			add(file);
		}
	}

	public void remove(File file) {

		MultiStateConfSpace confSpace = TestCoffee.affinity_2RL0_7mut();
		var state = confSpace.getState("complex");

		var db = new FixedDB(file, 2*1024*1024);
		var index = new NodeIndex(db, "index", state, null);

		// add a bunch of random nodes, fill up all the space
		Random rand = new Random(12345);
		for (int i=0; i<50_000; i++) {
			index.add(new NodeIndex.Node(
				state.index,
				Conf.make(state.confSpace),
				new BigExp(rand.nextDouble(), rand.nextInt())
			));
		}

		// poll some nodes
		long size = index.size();
		long numPolled = 10;
		for (int i=0; i<numPolled; i++) {
			index.remove(index.highestScore());
		}

		assertThat(index.size(), is(size - numPolled));
	}

	@Test
	public void remove_mem() {
		remove(null);
	}

	@Test
	public void remove_file() {
		try (var file = new TestBase.TempFile("node.index")) {
			remove(file);
		}
	}
}
