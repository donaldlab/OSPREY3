package edu.duke.cs.osprey.confspace;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.Test;

import java.util.*;


@SuppressWarnings("unchecked")
public class TestTupleTree {

	@Test
	public void single() {

		TupleTree<Double> tree = new TupleTree<>();
		tree.put(new RCTuple(0, 0), 5.0);

		assertThat(getEntries(tree, new int[] { 0 }), containsInAnyOrder(
			entry(new RCTuple(0, 0), 5.0)
		));

		assertThat(getEntries(tree, new int[] { -1 }), is(empty()));
		assertThat(getEntries(tree, new int[] { 1 }), is(empty()));
	}

	@Test
	public void singles() {

		TupleTree<Double> tree = new TupleTree<>();
		tree.put(new RCTuple(0, 0), 5.0);
		tree.put(new RCTuple(0, 2), 6.0);
		tree.put(new RCTuple(1, 0), 7.0);
		tree.put(new RCTuple(1, 3), 8.0);

		assertThat(getEntries(tree, new int[] { -1, -1 }), is(empty()));

		assertThat(getEntries(tree, new int[] { 0, 0 }), containsInAnyOrder(
			entry(new RCTuple(0, 0), 5.0),
			entry(new RCTuple(1, 0), 7.0)
		));

		assertThat(getEntries(tree, new int[] { 1, 1 }), is(empty()));

		assertThat(getEntries(tree, new int[] { 1, 0 }), containsInAnyOrder(
			entry(new RCTuple(1, 0), 7.0)
		));

		assertThat(getEntries(tree, new int[] { 0, 1 }), containsInAnyOrder(
			entry(new RCTuple(0, 0), 5.0)
		));

		assertThat(getEntries(tree, new int[] { 2, 3 }), containsInAnyOrder(
			entry(new RCTuple(0, 2), 6.0),
			entry(new RCTuple(1, 3), 8.0)
		));

		assertThat(getEntries(tree, new int[] { -1, 3 }), containsInAnyOrder(
			entry(new RCTuple(1, 3), 8.0)
		));

		assertThat(getEntries(tree, new int[] { 2, -1 }), containsInAnyOrder(
			entry(new RCTuple(0, 2), 6.0)
		));
	}

	@Test
	public void singlesWithRC1() {

		TupleTree<Double> tree = new TupleTree<>();
		tree.put(new RCTuple(0, 0), 5.0);
		tree.put(new RCTuple(0, 2), 6.0);
		tree.put(new RCTuple(1, 0), 7.0);
		tree.put(new RCTuple(1, 3), 8.0);

		assertThat(getEntriesWithPos1(tree, 0, new int[] { 0, 0 }), containsInAnyOrder(
			entry(new RCTuple(0, 0), 5.0)
		));
		assertThat(getEntriesWithPos1(tree, 1, new int[] { 0, 0 }), containsInAnyOrder(
			entry(new RCTuple(1, 0), 7.0)
		));

		assertThat(getEntriesWithPos1(tree, 0, new int[] { 1, 1 }), is(empty()));
		assertThat(getEntriesWithPos1(tree, 0, new int[] { 1, 1 }), is(empty()));
		assertThat(getEntriesWithPos1(tree, 1, new int[] { 1, 1 }), is(empty()));
		assertThat(getEntriesWithPos1(tree, 1, new int[] { 1, 1 }), is(empty()));

		assertThat(getEntriesWithPos1(tree, 1, new int[] { 1, 0 }), containsInAnyOrder(
			entry(new RCTuple(1, 0), 7.0)
		));
		assertThat(getEntriesWithPos1(tree, 0, new int[] { 1, 0 }), is(empty()));

		assertThat(getEntriesWithPos1(tree, 0, new int[] { 0, 1 }), containsInAnyOrder(
			entry(new RCTuple(0, 0), 5.0)
		));
		assertThat(getEntriesWithPos1(tree, 1, new int[] { 0, 1 }), is(empty()));

		assertThat(getEntriesWithPos1(tree, 0, new int[] { 2, 3 }), containsInAnyOrder(
			entry(new RCTuple(0, 2), 6.0)
		));
		assertThat(getEntriesWithPos1(tree, 1, new int[] { 2, 3 }), containsInAnyOrder(
			entry(new RCTuple(1, 3), 8.0)
		));
	}

	@Test
	public void pair() {

		TupleTree<Double> tree = new TupleTree<>();
		tree.put(new RCTuple(0, 0, 1, 0), 5.0);

		assertThat(getEntries(tree, new int[] { 0 }), is(empty()));

		assertThat(getEntries(tree, new int[] { 0, 0 }), containsInAnyOrder(
			entry(new RCTuple(0, 0, 1, 0), 5.0)
		));

		assertThat(getEntries(tree, new int[] { 0, 1 }), is(empty()));
		assertThat(getEntries(tree, new int[] { 1, 0 }), is(empty()));
		assertThat(getEntries(tree, new int[] { 1, 1 }), is(empty()));
	}

	@Test
	public void pairs() {

		TupleTree<Double> tree = new TupleTree<>();
		tree.put(new RCTuple(0, 0, 1, 0), 5.0);
		tree.put(new RCTuple(0, 0, 1, 1), 6.0);
		tree.put(new RCTuple(0, 1, 1, 0), 7.0);
		tree.put(new RCTuple(0, 0, 2, 0), 8.0);
		tree.put(new RCTuple(0, 0, 2, 1), 9.0);
		tree.put(new RCTuple(1, 1, 2, 0), 10.0);

		assertThat(getEntries(tree, new int[] { 0 }), is(empty()));
		assertThat(getEntries(tree, new int[] { 0, 2 }), is(empty()));

		assertThat(getEntries(tree, new int[] { 0, 0 }), containsInAnyOrder(
			entry(new RCTuple(0, 0, 1, 0), 5.0)
		));

		assertThat(getEntries(tree, new int[] { 0, 0, 0 }), containsInAnyOrder(
			entry(new RCTuple(0, 0, 1, 0), 5.0),
			entry(new RCTuple(0, 0, 2, 0), 8.0)
		));

		assertThat(getEntries(tree, new int[] { 0, 0, 1 }), containsInAnyOrder(
			entry(new RCTuple(0, 0, 1, 0), 5.0),
			entry(new RCTuple(0, 0, 2, 1), 9.0)
		));

		assertThat(getEntries(tree, new int[] { 0, 1, 0 }), containsInAnyOrder(
			entry(new RCTuple(0, 0, 1, 1), 6.0),
			entry(new RCTuple(0, 0, 2, 0), 8.0),
			entry(new RCTuple(1, 1, 2, 0), 10.0)
		));

		assertThat(getEntries(tree, new int[] { 0, 1, 1 }), containsInAnyOrder(
			entry(new RCTuple(0, 0, 1, 1), 6.0),
			entry(new RCTuple(0, 0, 2, 1), 9.0)
		));

		assertThat(getEntries(tree, new int[] { 1, 0, 0 }), containsInAnyOrder(
			entry(new RCTuple(0, 1, 1, 0), 7.0)
		));

		assertThat(getEntries(tree, new int[] { 1, 0, 1 }), containsInAnyOrder(
			entry(new RCTuple(0, 1, 1, 0), 7.0)
		));

		assertThat(getEntries(tree, new int[] { 1, 1, 0 }), containsInAnyOrder(
			entry(new RCTuple(1, 1, 2, 0), 10.0)
		));

		assertThat(getEntries(tree, new int[] { 1, 1, 1 }), is(empty()));

		assertThat(getEntries(tree, new int[] { 2, 0, 0 }), is(empty()));
		assertThat(getEntries(tree, new int[] { 0, 2, 0 }), containsInAnyOrder(
			entry(new RCTuple(0, 0, 2, 0), 8.0)
		));
		assertThat(getEntries(tree, new int[] { 0, 0, 2 }), containsInAnyOrder(
			entry(new RCTuple(0, 0, 1, 0), 5.0)
		));
	}

	// TODO: pairs with Pos1 and (Pos1 and Pos2)


	private static <T> Map.Entry<RCTuple,T> entry(RCTuple tuple, T value) {
		return new AbstractMap.SimpleEntry<>(tuple, value);
	}

	private static <T> List<Map.Entry<RCTuple,T>> getEntries(TupleTree<T> tree, int[] conf) {

		List<Map.Entry<RCTuple,T>> entries = new ArrayList<>();
		Set<RCTuple> tuples = new HashSet<>();

		tree.forEachIn(conf, (tuple, value) -> {

			// make sure we see each tuple only once
			assertThat(tuples.contains(tuple), is(false));
			tuples.add(tuple);

			entries.add(entry(tuple, value));
		});

		return entries;
	}

	private static <T> List<Map.Entry<RCTuple,T>> getEntriesWithPos1(TupleTree<T> tree, int pos1, int[] conf) {

		List<Map.Entry<RCTuple,T>> entries = new ArrayList<>();
		Set<RCTuple> tuples = new HashSet<>();

		tree.forEachIn(conf, pos1, (tuple, value) -> {

			// make sure we see each tuple only once
			assertThat(tuples.contains(tuple), is(false));
			tuples.add(tuple);

			entries.add(entry(tuple, value));
		});

		return entries;
	}
}
