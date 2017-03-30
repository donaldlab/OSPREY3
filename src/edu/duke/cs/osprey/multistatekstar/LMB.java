package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.util.StringTokenizer;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * generates a linear combination of partition function values over all states
 * constTerm + sum_s ( coeff_s * partition_function_value_for_state_s)
 * coeffs: >0 is upper bound, <0 is lower bound
 * interpretation: constraint LMBs are accepted iff eval() is <=0
 */
public class LMB {

	BigDecimal[] coeffs;
	BigDecimal constTerm;
	boolean correctConstTerm;

	public LMB(String s, int numStates) {
		StringTokenizer st = new StringTokenizer(s);
		if(st.countTokens()!=numStates+1){
			throw new RuntimeException("ERROR: there are "+numStates+" states but LinFunc "
					+ "specified with "+st.countTokens()+" coefficients: "+s);
		}

		coeffs = new BigDecimal[numStates];
		for(int state=0; state<numStates; state++)
			coeffs[state] = new BigDecimal( st.nextToken() );

		String s1 = st.nextToken().replaceAll("\\s+","").toLowerCase();   		
		if(s1.endsWith("xwt")) {
			constTerm = new BigDecimal(s1.split("x")[0]);
			correctConstTerm = true;
		}
		else {
			constTerm = new BigDecimal(s1);
			correctConstTerm = false;
		}
	}

	BigDecimal[] getCoeffs() {
		return coeffs;
	}

	boolean correctConstTerm() {
		return correctConstTerm;
	}

	void correctConstTerm(BigDecimal correction) {
		constTerm = constTerm.multiply(correction);
	}

	BigDecimal getConstTerm() {
		return constTerm;
	}

	BigDecimal eval(BigDecimal[] stateVals) {
		if(stateVals.length != coeffs.length)
			throw new RuntimeException("ERROR: Wrong number of state values");

		BigDecimal ans = BigDecimal.ZERO;

		for(int c=0; c<coeffs.length; c++){
			if(coeffs[c].compareTo(BigDecimal.ZERO) != 0)
				ans = ans.add(coeffs[c].multiply(stateVals[c]));
		}

		ans = ans.add(constTerm);
		return ans;
	}

}
