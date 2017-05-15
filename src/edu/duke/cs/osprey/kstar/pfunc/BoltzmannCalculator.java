package edu.duke.cs.osprey.kstar.pfunc;

import java.io.Serializable;
import java.math.BigDecimal;

import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.tools.ExpFunction;

@SuppressWarnings("serial")
public class BoltzmannCalculator implements Serializable {
	
	private double constRT = PoissonBoltzmannEnergy.constRT;
	private ExpFunction e = new ExpFunction();
	
	public BigDecimal calc(double energy) {
		return e.exp(-energy/constRT);
	}
}
