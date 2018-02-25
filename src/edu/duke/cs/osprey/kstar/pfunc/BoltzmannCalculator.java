package edu.duke.cs.osprey.kstar.pfunc;

import java.math.BigDecimal;
import java.math.MathContext;

import ch.obermuhlner.math.big.BigDecimalMath;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.tools.ExpFunction;

public class BoltzmannCalculator {
	
	public static double constRT = PoissonBoltzmannEnergy.constRT;

	private MathContext mathContext;

	public BoltzmannCalculator(MathContext mathContext) {
		this.mathContext = mathContext;
	}
	
	public BigDecimal calc(double energy) {
		return BigDecimalMath.exp(BigDecimal.valueOf(-energy/constRT), mathContext);
	}
}
