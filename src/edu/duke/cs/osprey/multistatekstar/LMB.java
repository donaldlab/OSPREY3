/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

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
