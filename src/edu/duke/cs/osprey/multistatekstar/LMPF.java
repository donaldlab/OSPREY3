package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;

public class LMPF {

	double[] coeffs;
    double constTerm;
	
	public LMPF(String s, int numStates) {
	}
	
	/**
	 * Default constructor that initializes coeffs and constant term 
	 * with default values.
	 * Default coeffs: 
	 * UbConstr00 -1 0 0 1e-8
	 * UbConstr01 0 -1 0 1e-8
	 */
	public LMPF() {
		
	}
	
    double[] getCoeffs() {
        return coeffs;
    }

    double getConstTerm() {
        return constTerm;
    }
    
    BigDecimal eval(BigDecimal[] statePFVals) {
    	return BigDecimal.ZERO;
    }
	
}
