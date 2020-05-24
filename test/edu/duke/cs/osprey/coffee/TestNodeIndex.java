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

		// track the number of evicted nodes
		AtomicLong numEvicted = new AtomicLong(0);
		Consumer<NodeIndex.StateConf> evictionListener = conf -> numEvicted.incrementAndGet();

		var index = new NodeIndex(file, 2*1024*1024, confSpace, evictionListener);

		// add a bunch of random nodes from the complex state
		var state = confSpace.getState("complex");
		Random rand = new Random(12345);
		for (int i=0; i<100_000; i++) {
			NodeIndex.StateConf sconf = new NodeIndex.StateConf(state.index, Conf.make(state.confSpace));
			BigExp zSumUpper = new BigExp(rand.nextDouble(), rand.nextInt());
			index.add(zSumUpper, sconf);
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

		var index = new NodeIndex(file, 2*1024*1024, confSpace, null);

		// add a bunch of random nodes from the complex state, fill up all the space
		var state = confSpace.getState("complex");
		Random rand = new Random(12345);
		for (int i=0; i<50_000; i++) {
			NodeIndex.StateConf sconf = new NodeIndex.StateConf(state.index, Conf.make(state.confSpace));
			BigExp zSumUpper = new BigExp(rand.nextDouble(), rand.nextInt());
			index.add(zSumUpper, sconf);
		}

		// poll some nodes
		long size = index.size();
		long numPolled = 10;
		for (int i=0; i<numPolled; i++) {
			index.remove(index.highestKey());
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
