package edu.duke.cs.osprey.tools;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class MathTools {
	
	public static int divUp(int num, int denom) {
		return (num + denom - 1)/denom;
	}
	
	public static int roundUpToMultiple(int val, int base) {
		int mod = val % base;
		if (mod == 0) {
			return val;
		}
		return val + base - mod;
	}

	public static <T> List<List<T>> powerset(List<T> list) {
		return powersetUpTo(list, list.size());
	}

	public static <T> List<List<T>> powersetUpTo(List<T> list, int size) {

		// adapted from powerset algorithm:
		// http://rosettacode.org/wiki/Power_set#Iterative

		List<List<T>> powerset = new ArrayList<>();

		// start with the empty set
		powerset.add(new ArrayList<>());

		for (T item : list) {

			for (List<T> subset : new ArrayList<>(powerset)) {

				// expand subset if within size
				if (subset.size() < size) {
					List<T> newSubset = new ArrayList<>(subset);
					newSubset.add(item);
					powerset.add(newSubset);
				}
			}
		}
		return powerset;
	}

	public static <T> Iterable<List<T>> cartesianProduct(List<List<T>> lists) {

		// stack overflow is awesome for the lazy programmer:
		// https://stackoverflow.com/questions/714108/cartesian-product-of-arbitrary-sets-in-java
		// I mean, sure I could figure out how to do this stuff from first principles, but why would I want to?

		return () -> {

			return new Iterator<List<T>>() {

				private final int n;
				private boolean hasNext;
				private final int[] lengths;
				private final int[] indices;

				{
					n = lists.size();
					hasNext = n > 0;

					lengths = new int[n];
					for (int i=0; i<n; i++) {
						lengths[i] = lists.get(i).size();
						if (lengths[i] <= 0) {
							hasNext = false;
						}
					}

					indices = new int[n];
					Arrays.fill(indices, 0);
				}

				public boolean hasNext() {
					return hasNext;
				}

				public List<T> next() {

					// build the output
					List<T> result = new ArrayList<>();
					for (int i=0; i<n; i++) {
						result.add(lists.get(i).get(indices[i]));
					}

					// advance the indices
					for (int i = n - 1; i >= 0; i--) {
						if (indices[i] == lengths[i] - 1) {
							indices[i] = 0;
							if (i == 0) {
								hasNext = false;
							}
						} else {
							indices[i]++;
							break;
						}
					}

					return result;
				}
			};
		};
	}
}
