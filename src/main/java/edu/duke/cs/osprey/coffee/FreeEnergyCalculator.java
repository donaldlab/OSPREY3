package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;


public class FreeEnergyCalculator {

	public final MathContext mathContext;
	public final BoltzmannCalculator bcalc;

	public FreeEnergyCalculator() {

		// only computing ln, so only need enough precision to fill a double
		mathContext = new MathContext(16, RoundingMode.HALF_UP);
		bcalc = new BoltzmannCalculator(mathContext);
	}

	public double calc(BigDecimal z) {
		return bcalc.freeEnergyPrecise(z);
	}

	public MathTools.DoubleBounds calc(MathTools.BigDecimalBounds z) {
		return bcalc.freeEnergyPrecise(z);
	}

	public double calc(BigExp z) {
		return bcalc.freeEnergyPrecise(z);
	}
}
