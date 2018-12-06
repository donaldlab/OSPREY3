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

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.tools.HashCalculator;

import java.util.*;

public class Conf {

	public static final int Unassigned = -1;

	public static int[] make(SimpleConfSpace confSpace) {
		int[] conf = new int[confSpace.positions.size()];
		Arrays.fill(conf, Unassigned);
		return conf;
	}

	public static int[] make(SimpleConfSpace confSpace, RCTuple frag) {

		// start with an unassigned conf
		int[] conf = make(confSpace);

		// assign the fragment to the conf
		for (int i=0; i<frag.size(); i++) {
			conf[frag.pos.get(i)] = frag.RCs.get(i);
		}

		return conf;
	}

	public static int[] make(ConfIndex index) {
		int[] conf = new int[index.numPos];
		Conf.clear(conf);
		for (int i=0; i<index.numDefined; i++) {
			conf[index.definedPos[i]] = index.definedRCs[i];
		}
		return conf;
	}

	public static void index(int[] conf, ConfIndex index) {
		index.numDefined = 0;
		index.numUndefined = 0;
		for (int pos=0; pos<conf.length; pos++) {
			int rc = conf[pos];
			if (rc == Unassigned) {
				index.undefinedPos[index.numUndefined] = pos;
				index.numUndefined++;
			} else {
				index.definedPos[index.numDefined] = pos;
				index.definedRCs[index.numDefined] = conf[pos];
				index.numDefined++;
			}
		}
		index.node = null;
	}

	public static ConfIndex index(int[] conf) {
		ConfIndex index = new ConfIndex(conf.length);
		index(conf, index);
		return index;
	}

	public static void clear(int[] conf) {
		Arrays.fill(conf, Unassigned);
	}

	public static int hashCode(int[] conf) {
		return Arrays.hashCode(conf);
	}

	public static boolean equals(int[] a, int[] b) {
		return Arrays.equals(a, b);
	}

	public static String toString(int[] conf) {
		return Arrays.toString(conf);
	}

	public static boolean isCompletelyAssigned(int[] conf) {
		for (int i : conf) {
			if (i == Unassigned) {
				return false;
			}
		}
		return true;
	}

	public static List<RCTuple> getPairs(int[] conf) {

		List<RCTuple> pairs = new ArrayList<>();

		for (int pos1=0; pos1<conf.length; pos1++) {

			int rc1 = conf[pos1];
			if (rc1 == Unassigned) {
				continue;
			}

			for (int pos2=0; pos2<pos1; pos2++) {

				int rc2 = conf[pos2];
				if (rc2 == Unassigned) {
					continue;
				}

				// NOTE: reverse order so the tuple is sorted, since pos2 < pos1
				pairs.add(new RCTuple(pos2, rc2, pos1, rc1));
			}
		}

		return pairs;
	}

	public static boolean containsTuple(int[] conf, RCTuple tuple) {
		for (int i=0; i<tuple.size(); i++) {
			int pos = tuple.pos.get(i);
			int rc = tuple.RCs.get(i);
			if (conf[pos] != rc) {
				return false;
			}
		}
		return true;
	}

	public static boolean contains(Collection<int[]> confs, int [] conf) {
		for (int[] c : confs) {
			if (equals(c, conf)) {
				return true;
			}
		}
		return false;
	}


	// sigh, this is the only way I know to attach custom hashers and equality checkers to a hash set
	private static class Wrapper {

		final int[] assignments;

		public Wrapper(int[] assignments) {
			this.assignments = assignments;
		}

		@Override
		public int hashCode() {
			return Conf.hashCode(assignments);
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof Wrapper && equals((Wrapper)other);
		}

		public boolean equals(Wrapper other) {
			return Conf.equals(this.assignments, other.assignments);
		}

		@Override
		public String toString() {
			return Conf.toString(assignments);
		}
	}

	private static int[] checkKey(Object o) {
		if (o == null) {
			throw new IllegalArgumentException("key can't be null");
		}
		if (o instanceof int[]) {
			return (int[])o;
		}
		throw new IllegalArgumentException("key is a " + o.getClass().getName() + " instead of an int[]");
	}

	public static class Set extends AbstractSet<int[]> {

		private java.util.Set<Wrapper> set;

		public Set() {
			this(new HashSet<Wrapper>());
		}

		public Set(java.util.Set<Wrapper> set) {
			this.set = set;
		}

		public Set(Collection<int[]> source) {
			this();
			addAll(source);
		}

		public Set(int[] ... confs) {
			this();
			for (int[] conf : confs) {
				add(conf);
			}
		}

		@Override
		public Iterator<int[]> iterator() {
			return new Iterator<int[]>() {

				Iterator<Wrapper> iter = set.iterator();

				@Override
				public boolean hasNext() {
					return iter.hasNext();
				}

				@Override
				public int[] next() {
					return iter.next().assignments;
				}

				@Override
				public void remove() {
					iter.remove();
				}
			};
		}

		@Override
		public int size() {
			return set.size();
		}

		@Override
		public boolean isEmpty() {
			return set.isEmpty();
		}

		@Override
		public boolean contains(Object o) {
			return contains(checkKey(o));
		}

		public boolean contains(int[] conf) {
			return set.contains(new Wrapper(conf));
		}

		@Override
		public boolean add(int[] conf) {
			return set.add(new Wrapper(conf));
		}

		@Override
		public boolean remove(Object o) {
			return remove(checkKey(o));
		}

		public boolean remove(int[] conf) {
			return set.remove(new Wrapper(conf));
		}

		@Override
		public void clear() {
			set.clear();
		}

		@Override
		public boolean equals(Object o) {
			return o instanceof Set && equals((Set)o);
		}

		public boolean equals(Set other) {
			return this.set.equals(other.set);
		}

		@Override
		public int hashCode() {
			return set.hashCode();
		}

		@Override
		public String toString() {
			return set.toString();
		}
	}

	public static class Map<V> extends AbstractMap<int[],V> {

		private final HashMap<Wrapper,V> map = new HashMap<>();

		@Override
		public int size() {
			return map.size();
		}

		@Override
		public boolean isEmpty() {
			return map.isEmpty();
		}

		@Override
		public boolean containsKey(Object key) {
			return key instanceof int[] && containsKey((int[])key);
		}

		public boolean containsKey(int[] key) {
			return map.containsKey(new Wrapper(key));
		}

		@Override
		public boolean containsValue(Object value) {
			return map.containsValue(value);
		}

		@Override
		public V get(Object key) {
			return get(checkKey(key));
		}

		public V get(int[] key) {
			return map.get(new Wrapper(key));
		}

		@Override
		public V put(int[] key, V value) {
			return map.put(new Wrapper(key), value);
		}

		@Override
		public V remove(Object key) {
			return remove(checkKey(key));
		}

		public V remove(int[] key) {
			return map.remove(new Wrapper(key));
		}

		@Override
		public void clear() {
			map.clear();
		}

		@Override
		public java.util.Set<int[]> keySet() {
			return new Conf.Set(map.keySet());
		}

		@Override
		public Collection<V> values() {
			return map.values();
		}

		private static class EntryWrapper<V> implements Entry<Wrapper,V> {

			private Entry<int[],V> entry;
			private Wrapper wrapper;

			public EntryWrapper(Entry<int[],V> entry) {
				this.entry = entry;
				this.wrapper = new Wrapper(entry.getKey());
			}

			@Override
			public Wrapper getKey() {
				return wrapper;
			}

			@Override
			public V getValue() {
				return entry.getValue();
			}

			@Override
			public V setValue(V value) {
				return entry.setValue(value);
			}

			@Override
			public int hashCode() {
				return HashCalculator.combineHashes(
					wrapper.hashCode(),
					Objects.hashCode(entry.getValue())
				);
			}

			@Override
			@SuppressWarnings("unchecked")
			public boolean equals(Object o) {
				return o instanceof Entry && equals((Entry)o);
			}

			public boolean equals(Entry<int[],V> other) {
				return Conf.equals(this.wrapper.assignments, other.getKey())
					&& Objects.equals(this.getValue(), other.getValue());
			}

			@Override
			public String toString() {
				return String.format("%s=%s", getKey(), getValue());
			}
		}

		private static class EntryUnwrapper<V> implements Entry<int[],V> {

			private Entry<Wrapper,V> entry;

			public EntryUnwrapper(Entry<Wrapper,V> entry) {
				this.entry = entry;
			}

			@Override
			public int[] getKey() {
				return entry.getKey().assignments;
			}

			@Override
			public V getValue() {
				return entry.getValue();
			}

			@Override
			public V setValue(V value) {
				return entry.setValue(value);
			}

			@Override
			public int hashCode() {
				return entry.hashCode();
			}

			@Override
			@SuppressWarnings("unchecked")
			public boolean equals(Object o) {
				return o instanceof Entry && equals((Entry)o);
			}

			public boolean equals(Entry<int[],V> other) {
				return Conf.equals(this.getKey(), other.getKey())
					&& Objects.equals(this.getValue(), other.getValue());
			}

			@Override
			public String toString() {
				return String.format("%s=%s", entry.getKey(), getValue());
			}
		}

		@Override
		public java.util.Set<Entry<int[],V>> entrySet() {

			return new AbstractSet<Entry<int[],V>>() {

				private final java.util.Set<Entry<Wrapper,V>> set = map.entrySet();

				@Override
				public Iterator<Entry<int[],V>> iterator() {
					return new Iterator<Entry<int[],V>>() {

						Iterator<Entry<Wrapper,V>> iter = set.iterator();

						@Override
						public boolean hasNext() {
							return iter.hasNext();
						}

						@Override
						public Entry<int[],V> next() {
							return new EntryUnwrapper<>(iter.next());
						}

						@Override
						public void remove() {
							iter.remove();
						}
					};
				}

				@Override
				public int size() {
					return set.size();
				}

				@Override
				public boolean isEmpty() {
					return set.isEmpty();
				}

				@Override
				@SuppressWarnings("unchecked")
				public boolean contains(Object o) {
					return o instanceof Entry && contains((Entry)o);
				}

				public boolean contains(Entry<int[],V> entry) {
					return set.contains(new EntryWrapper<>(entry));
				}

				@Override
				public boolean add(Entry<int[],V> entry) {
					return set.add(new EntryWrapper<>(entry));
				}

				@Override
				@SuppressWarnings("unchecked")
				public boolean remove(Object o) {
					return o instanceof Entry && remove((Entry)o);
				}

				public boolean remove(Entry<int[],V> entry) {
					return set.remove(new EntryWrapper<>(entry));
				}

				@Override
				public void clear() {
					set.clear();
				}

				@Override
				public boolean equals(Object o) {
					return o instanceof Set && equals((Set)o);
				}

				public boolean equals(Set other) {
					return this.set.equals(other.set);
				}

				@Override
				public int hashCode() {
					return set.hashCode();
				}

				@Override
				public String toString() {
					return set.toString();
				}
			};
		}

		@Override
		public boolean equals(Object o) {
			return o instanceof Map && equals((Map)o);
		}

		public boolean equals(Map other) {
			return this.map.equals(other.map);
		}

		@Override
		public int hashCode() {
			return map.hashCode();
		}

		@Override
		public String toString() {
			return map.toString();
		}
	}
}
