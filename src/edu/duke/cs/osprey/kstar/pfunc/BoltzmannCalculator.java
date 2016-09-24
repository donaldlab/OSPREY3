package edu.duke.cs.osprey.kstar.pfunc;

import java.math.BigDecimal;

import edu.duke.cs.osprey.tools.ExpFunction;

public class BoltzmannCalculator {
	
	public static final double RT = 1.9891/1000 * 298.15;
	
	private ExpFunction e = new ExpFunction();
	
	public BigDecimal calc(double energy) {
		return e.exp(-energy/RT);
	}
}
