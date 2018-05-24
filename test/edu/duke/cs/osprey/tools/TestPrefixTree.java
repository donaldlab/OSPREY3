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

package edu.duke.cs.osprey.tools;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.Test;

import java.util.Arrays;

public class TestPrefixTree {

	@Test
	public void empty() {

		PrefixTreeSet<Integer> tree = new PrefixTreeSet<>();

		assertThat(tree.contains(Arrays.asList()), is(false));

		assertThat(tree.contains(Arrays.asList(1)), is(false));
		assertThat(tree.contains(Arrays.asList(1, 2, 3)), is(false));
	}

	@Test
	public void one() {

		PrefixTreeSet<Integer> tree = new PrefixTreeSet<>();

		tree.add(Arrays.asList(1));

		assertThat(tree.contains(Arrays.asList()), is(false));

		assertThat(tree.contains(Arrays.asList(1)), is(true));

		assertThat(tree.contains(Arrays.asList(2)), is(false));
		assertThat(tree.contains(Arrays.asList(1, 2, 3)), is(false));
	}

	@Test
	public void oneTwo() {

		PrefixTreeSet<Integer> tree = new PrefixTreeSet<>();

		tree.add(Arrays.asList(1, 2));

		assertThat(tree.contains(Arrays.asList()), is(false));
		assertThat(tree.contains(Arrays.asList(1)), is(false));

		assertThat(tree.contains(Arrays.asList(1, 2)), is(true));

		assertThat(tree.contains(Arrays.asList(2)), is(false));
		assertThat(tree.contains(Arrays.asList(1, 2, 3)), is(false));
	}

	@Test
	public void multipleSequences() {

		PrefixTreeSet<Integer> tree = new PrefixTreeSet<>();

		tree.add(Arrays.asList(1, 2));
		tree.add(Arrays.asList(2, 3, 4));
		tree.add(Arrays.asList(9));
		tree.add(Arrays.asList(7, 2));

		assertThat(tree.contains(Arrays.asList(1)), is(false));
		assertThat(tree.contains(Arrays.asList(1, 2)), is(true));

		assertThat(tree.contains(Arrays.asList(2)), is(false));
		assertThat(tree.contains(Arrays.asList(2, 3)), is(false));
		assertThat(tree.contains(Arrays.asList(2, 3, 4)), is(true));

		assertThat(tree.contains(Arrays.asList(3)), is(false));
		assertThat(tree.contains(Arrays.asList(3, 1)), is(false));
		assertThat(tree.contains(Arrays.asList(3, 2)), is(false));
		assertThat(tree.contains(Arrays.asList(3, 2, 3, 4)), is(false));
		assertThat(tree.contains(Arrays.asList(3, 9)), is(false));

		assertThat(tree.contains(Arrays.asList(9)), is(true));

		assertThat(tree.contains(Arrays.asList(7)), is(false));
		assertThat(tree.contains(Arrays.asList(7, 2)), is(true));
	}
}
