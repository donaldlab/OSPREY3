package edu.duke.cs.osprey.energy;

public class ScaledEnergyFunction implements EnergyFunction {
	
	private static final long serialVersionUID = -1352051347168314110L;
	
	private EnergyFunction efunc;
	private double c;
	
	public ScaledEnergyFunction(EnergyFunction efunc, double c) {
		this.efunc = efunc;
		this.c = c;
	}

	@Override
	public double getEnergy() {
		return efunc.getEnergy()*c;
	}
}
