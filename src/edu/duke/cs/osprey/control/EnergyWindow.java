package edu.duke.cs.osprey.control;

public class EnergyWindow {
	
	private double min;
	private double size;
	
	public EnergyWindow(double energy, double size) {
		min = energy;
		this.size = size;
	}
	
	public boolean update(double energy) {
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
}
