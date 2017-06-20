package edu.duke.cs.osprey.gmec;

public class EnergyRange {
	
	private double min;
	private double size;
	
	public EnergyRange(double energy, double size) {
		min = energy;
		this.size = size;
	}
	
	public boolean updateMin(double energy) {
		if (energy >= min) {
			return false;
		}
		min = energy;
		return true;
	}
	
	public double getMin() {
		return min;
	}
	
	public double getMax() {
		return min + size;
	}

	public boolean contains(double energy) {
		return energy >= min && energy <= getMax();
	}
	
	public boolean containsOrBelow(double energy) {
		return energy <= getMax();
	}
}
