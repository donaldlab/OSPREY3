package edu.duke.cs.osprey.restypes;

import java.io.Serializable;
import java.util.HashMap;

public class ResidueEntropies implements Serializable {

	private HashMap<String,Double> values = new HashMap<>();

	public void set(String resType, double val) {
		values.put(resType, val);
	}

	public double get(String resType) {

		// return defined entropy, if any
		Double val = values.get(normalize(resType));
		if (val != null) {
			return val;
		}

		// or zero entropy, by default
		return 0;
	}

	public void clear() {
		values.clear();
	}

	public void setAll(ResidueEntropies other) {
		values.putAll(other.values);
	}

	private String normalize(String resType) {
		return resType.toUpperCase();
	}
}
