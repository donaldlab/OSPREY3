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

package edu.duke.cs.osprey.confspace;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.Test;

import java.util.*;


public class TestConf {

	@Test
	public void emptySet() {

		Conf.Set set = new Conf.Set();

		assertThat(set.isEmpty(), is(true));
		assertThat(set.size(), is(0));
		assertThat(set.iterator().hasNext(), is(false));
		assertThat(set.contains(new int[] { 0, 0, 0 }), is(false));
		assertThat(set, is(new Conf.Set()));
	}

	@Test
	public void oneConfSet() {

		Conf.Set set = new Conf.Set();

		set.add(new int[] { 0, 0, 0 });

		assertThat(set.isEmpty(), is(false));
		assertThat(set.size(), is(1));
		assertThat(set, containsInAnyOrder(new int[] { 0, 0, 0 }));
		assertThat(set.contains(new int[] { 0, 0, 0 }), is(true));
		assertThat(set, is(new Conf.Set(new int[] { 0, 0, 0 })));
	}

	@Test
	public void threeConfsSet() {

		Conf.Set set = new Conf.Set();

		set.add(new int[] { 0, 0, 0 });
		set.add(new int[] { 1, 2, 3 });
		set.add(new int[] { 9, 8, 7 });

		assertThat(set.isEmpty(), is(false));
		assertThat(set.size(), is(3));
		assertThat(set, containsInAnyOrder(
			new int[] { 0, 0, 0 },
			new int[] { 1, 2, 3 },
			new int[] { 9, 8, 7 }
		));
		assertThat(set.contains(new int[] { 0, 0, 0 }), is(true));
		assertThat(set.contains(new int[] { 1, 2, 3 }), is(true));
		assertThat(set.contains(new int[] { 9, 8, 7 }), is(true));
		assertThat(set, is(new Conf.Set(
			new int[] { 0, 0, 0 },
			new int[] { 1, 2, 3 },
			new int[] { 9, 8, 7 }
		)));
	}

	@Test
	public void setRemove() {

		Conf.Set set = new Conf.Set();

		set.add(new int[] { 0, 0, 0 });
		set.add(new int[] { 1, 2, 3 });
		set.add(new int[] { 9, 8, 7 });

		set.remove(new int[] { 0, 0, 0 });

		assertThat(set.isEmpty(), is(false));
		assertThat(set.size(), is(2));
		assertThat(set, containsInAnyOrder(
			new int[] { 1, 2, 3 },
			new int[] { 9, 8, 7 }
		));
		assertThat(set.contains(new int[] { 0, 0, 0 }), is(false));
		assertThat(set.contains(new int[] { 1, 2, 3 }), is(true));
		assertThat(set.contains(new int[] { 9, 8, 7 }), is(true));
		assertThat(set, is(new Conf.Set(
			new int[] { 1, 2, 3 },
			new int[] { 9, 8, 7 }
		)));
	}

	@Test
	public void duplicatesSet() {

		Conf.Set set = new Conf.Set();

		set.add(new int[] { 0, 0, 0 });
		set.add(new int[] { 0, 0, 0 });
		set.add(new int[] { 0, 0, 0 });

		assertThat(set.isEmpty(), is(false));
		assertThat(set.size(), is(1));
		assertThat(set, containsInAnyOrder(new int[] { 0, 0, 0 }));
		assertThat(set.contains(new int[] { 0, 0, 0 }), is(true));
		assertThat(set, is(new Conf.Set(new int[] { 0, 0, 0 })));
	}

	@Test
	public void emptyMap() {

		Conf.Map<Integer> map = new Conf.Map<>();

		assertThat(map.isEmpty(), is(true));
		assertThat(map.size(), is(0));
		assertThat(map.containsKey(new int[] { 0, 0, 0 }), is(false));
		assertThat(map.containsValue(42), is(false));

		assertThat(map.keySet().isEmpty(), is(true));
		assertThat(map.keySet().size(), is(0));
		assertThat(map.keySet().contains(new int[] { 0, 0, 0 }), is(false));

		assertThat(map.values().isEmpty(), is(true));
		assertThat(map.values().size(), is(0));

		assertThat(map.entrySet().isEmpty(), is(true));
		assertThat(map.entrySet().size(), is(0));
	}

	public static class TestEntry<V> extends AbstractMap.SimpleEntry<int[],V> {

		public TestEntry(int[] conf, V val) {
			super(conf, val);
		}

		@Override
		public String toString() {
			return String.format("%s=%s", Conf.toString(getKey()), getValue());
		}
	}

	@Test
	@SuppressWarnings("unchecked")
	public void oneMap() {

		Conf.Map<Integer> map = new Conf.Map<>();

		map.put(new int[] { 0, 0, 0 }, 5);

		assertThat(map.isEmpty(), is(false));
		assertThat(map.size(), is(1));
		assertThat(map.containsKey(new int[] { 0, 0, 0 }), is(true));
		assertThat(map.containsValue(5), is(true));

		assertThat(map.keySet().isEmpty(), is(false));
		assertThat(map.keySet().size(), is(1));
		assertThat(map.keySet().contains(new int[] { 0, 0, 0 }), is(true));
		assertThat(map.keySet(), containsInAnyOrder(new int[] { 0, 0, 0 }));

		assertThat(map.values().isEmpty(), is(false));
		assertThat(map.values().size(), is(1));
		assertThat(map.values(), containsInAnyOrder(5));

		assertThat(map.entrySet().isEmpty(), is(false));
		assertThat(map.entrySet().size(), is(1));
		assertThat(map.entrySet(), containsInAnyOrder(
			new TestEntry<>(new int[] { 0, 0, 0 }, 5)
		));
	}

	@Test
	@SuppressWarnings("unchecked")
	public void threeMap() {

		Conf.Map<Integer> map = new Conf.Map<>();

		map.put(new int[] { 0, 0, 0 }, 5);
		map.put(new int[] { 1, 2, 3 }, 7);
		map.put(new int[] { 9, 8, 7 }, 9);

		assertThat(map.isEmpty(), is(false));
		assertThat(map.size(), is(3));
		assertThat(map.containsKey(new int[] { 0, 0, 0 }), is(true));
		assertThat(map.containsKey(new int[] { 1, 2, 3 }), is(true));
		assertThat(map.containsKey(new int[] { 9, 8, 7 }), is(true));
		assertThat(map.containsValue(5), is(true));
		assertThat(map.containsValue(7), is(true));
		assertThat(map.containsValue(9), is(true));

		assertThat(map.keySet().isEmpty(), is(false));
		assertThat(map.keySet().size(), is(3));
		assertThat(map.keySet().contains(new int[] { 0, 0, 0 }), is(true));
		assertThat(map.keySet().contains(new int[] { 1, 2, 3 }), is(true));
		assertThat(map.keySet().contains(new int[] { 9, 8, 7 }), is(true));
		assertThat(map.keySet(), containsInAnyOrder(
			new int[] { 0, 0, 0 },
			new int[] { 1, 2, 3 },
			new int[] { 9, 8, 7 }
		));

		assertThat(map.values().isEmpty(), is(false));
		assertThat(map.values().size(), is(3));
		assertThat(map.values(), containsInAnyOrder(5, 7, 9));

		assertThat(map.entrySet().isEmpty(), is(false));
		assertThat(map.entrySet().size(), is(3));
		assertThat(map.entrySet(), containsInAnyOrder(
			new TestEntry<>(new int[] { 0, 0, 0 }, 5),
			new TestEntry<>(new int[] { 1, 2, 3 }, 7),
			new TestEntry<>(new int[] { 9, 8, 7 }, 9)
		));
	}

	@Test
	@SuppressWarnings("unchecked")
	public void mapRemove() {

		Conf.Map<Integer> map = new Conf.Map<>();

		map.put(new int[] { 0, 0, 0 }, 5);
		map.put(new int[] { 1, 2, 3 }, 7);
		map.put(new int[] { 9, 8, 7 }, 9);

		map.remove(new int[] { 0, 0, 0 });

		assertThat(map.isEmpty(), is(false));
		assertThat(map.size(), is(2));
		assertThat(map.containsKey(new int[] { 1, 2, 3 }), is(true));
		assertThat(map.containsKey(new int[] { 9, 8, 7 }), is(true));
		assertThat(map.containsValue(7), is(true));
		assertThat(map.containsValue(9), is(true));

		assertThat(map.keySet().isEmpty(), is(false));
		assertThat(map.keySet().size(), is(2));
		assertThat(map.keySet().contains(new int[] { 1, 2, 3 }), is(true));
		assertThat(map.keySet().contains(new int[] { 9, 8, 7 }), is(true));
		assertThat(map.keySet(), containsInAnyOrder(
			new int[] { 1, 2, 3 },
			new int[] { 9, 8, 7 }
		));

		assertThat(map.values().isEmpty(), is(false));
		assertThat(map.values().size(), is(2));
		assertThat(map.values(), containsInAnyOrder(7, 9));

		assertThat(map.entrySet().isEmpty(), is(false));
		assertThat(map.entrySet().size(), is(2));
		assertThat(map.entrySet(), containsInAnyOrder(
			new TestEntry<>(new int[] { 1, 2, 3 }, 7),
			new TestEntry<>(new int[] { 9, 8, 7 }, 9)
		));
	}

	@Test
	@SuppressWarnings("unchecked")
	public void duplicatesMap() {

		Conf.Map<Integer> map = new Conf.Map<>();

		map.put(new int[] { 0, 0, 0 }, 5);
		map.put(new int[] { 0, 0, 0 }, 7);
		map.put(new int[] { 0, 0, 0 }, 9);

		assertThat(map.isEmpty(), is(false));
		assertThat(map.size(), is(1));
		assertThat(map.containsKey(new int[] { 0, 0, 0 }), is(true));
		assertThat(map.containsValue(5), is(false));
		assertThat(map.containsValue(7), is(false));
		assertThat(map.containsValue(9), is(true));

		assertThat(map.keySet().isEmpty(), is(false));
		assertThat(map.keySet().size(), is(1));
		assertThat(map.keySet().contains(new int[] { 0, 0, 0 }), is(true));
		assertThat(map.keySet(), containsInAnyOrder(new int[] { 0, 0, 0 }));

		assertThat(map.values().isEmpty(), is(false));
		assertThat(map.values().size(), is(1));
		assertThat(map.values(), containsInAnyOrder(9));

		assertThat(map.entrySet().isEmpty(), is(false));
		assertThat(map.entrySet().size(), is(1));
		assertThat(map.entrySet(), containsInAnyOrder(
			new TestEntry<>(new int[] { 0, 0, 0 }, 9)
		));
	}
}
