package edu.duke.cs.osprey.coffee.directors;

import java.util.*;
import java.util.function.Function;


/**
 * Track the K best things by lowest value, allowing ties.
 */
public class BestSet<T> {

	public final int capacity;
	public final Function<T,Double> scorer;

	private final SortedMap<Double,List<T>> values = new TreeMap<>();
	private int size = 0;

	public BestSet(int capacity, Function<T,Double> scorer) {
		this.capacity = capacity;
		this.scorer = scorer;
	}

	public int size() {
		return size;
	}

	public void clear() {
		values.clear();
		size = 0;
	}

	public void add(T value) {

		// add the values to the list in the right order, allowing ties
		values
			.computeIfAbsent(scorer.apply(value), key -> new ArrayList<>())
			.add(value);
		size += 1;

		// if we're over the capacity, get rid of the worst-scoring value(s)
		if (size > capacity) {
			double worstScore = values.lastKey();
			var list = values.get(worstScore);
			if (size - list.size() >= capacity) {
				values.remove(worstScore);
				size -= list.size();
			}
			// unless getting rid of the worst-scoring values puts us under the capacity, then don't do it
		}
	}

	public double minScore() {
		return values.firstKey();
	}
	public List<T> minValues() {
		return values.get(values.firstKey());
	}

	public double maxScore() {
		return values.lastKey();
	}
	public List<T> maxValues() {
		return values.get(values.lastKey());
	}

	public boolean contains(T value) {
		var list = values.get(scorer.apply(value));
		return list != null && list.contains(value);
	}

	public List<T> toList() {
		var list = new ArrayList<T>(size);
		values.forEach((s, v) -> list.addAll(v));
		return list;
	}

	public Iterable<Double> scores() {
		return values.keySet();
	}

	public List<T> values(double score) {
		return values.get(score);
	}
}
