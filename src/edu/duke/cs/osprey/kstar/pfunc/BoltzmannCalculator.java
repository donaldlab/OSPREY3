package edu.duke.cs.osprey.kstar.pfunc;

import java.math.BigDecimal;
import java.math.MathContext;

import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.tools.ExpFunction;

public class BoltzmannCalculator {
	
	public static double constRT = PoissonBoltzmannEnergy.constRT;

	public final ExpFunction e;

	public BoltzmannCalculator(MathContext mathContext) {
		e = new ExpFunction(mathContext);
	}
	
	public BigDecimal calc(double energy) {
		return e.exp(-energy/constRT);
	}
}
