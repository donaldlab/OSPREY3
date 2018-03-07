package edu.duke.cs.osprey.confspace;

import java.io.Serializable;
import java.util.*;

public class Conf {

	public static final int Unassigned = -1;

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
			this(new HashSet<>());
		}

		public Set(java.util.Set<Wrapper> set) {
			this.set = set;
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
			return contains(checkKey(0));
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
							Entry<Wrapper,V> entry = iter.next();
							Wrapper wrapper = entry.getKey();
							return new Entry<int[],V>() {

								@Override
								public int[] getKey() {
									return wrapper.assignments;
								}

								@Override
								public V getValue() {
									return entry.getValue();
								}

								@Override
								public V setValue(V value) {
									return entry.setValue(value);
								}
							};
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

				// TODO: NEXTTIME: finish implementing this interface

				@Override
				public boolean contains(Object o) {
					return o instanceof int[] && set.contains(new Wrapper((int[])o));
				}

				@Override
				public boolean add(int[] conf) {
					return set.add(new Wrapper(conf));
				}

				@Override
				public boolean remove(Object o) {
					return o instanceof int[] && set.remove(new Wrapper((int[])o));
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

				// TODO: until here
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
	}
}
