/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
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
