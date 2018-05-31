package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;

import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;

@SuppressWarnings("serial")
public class MSBoltzmannCalculator extends BoltzmannCalculator {

	private final static double MIN_ENERGY = -4000;
	
	@Override
	public BigDecimal calc(double energy) {
		// avoid overflow by capping min energy
		if(energy < MIN_ENERGY) energy = MIN_ENERGY;
		return super.calc(energy);
	}
	
}
