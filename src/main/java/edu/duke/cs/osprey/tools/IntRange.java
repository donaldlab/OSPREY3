package edu.duke.cs.osprey.tools;


public class IntRange {

	/** inclusive */
	public final int min;

	/** exclusive */
	public final int max;

	public IntRange(int min, int max) {
		this.min = min;
		this.max = max;
	}

	public boolean contains(int i) {
		return i >= min && i < max;
	}

	public int size() {
		return max - min;
	}
}
