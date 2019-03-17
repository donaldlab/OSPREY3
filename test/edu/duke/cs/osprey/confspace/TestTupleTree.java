/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

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
