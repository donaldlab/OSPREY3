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

package edu.duke.cs.osprey.kstar;

import java.math.BigDecimal;
import java.math.BigInteger;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;

@SuppressWarnings("serial")
public class KSQPrimeCalculator extends KSPStarCalculator {

	private BigInteger unPruned = BigInteger.ZERO;

	public KSQPrimeCalculator( PFAbstract pf ) {
		super(pf);

		unPruned = pf.getNumUnPruned();
		
		confSearch = getConfSearch();
	}

	
	private ConfSearch getConfSearch() {
		return pf.getConfTree(false);
	}
	

	public BigDecimal getQPrime( BigDecimal partialQLB ) {

		BigDecimal ans = BigDecimal.ZERO;

		synchronized( lock ) {
			ans = totalPF.subtract(partialQLB);

			BigDecimal uniformBound = new BigDecimal(unPruned.subtract(enumerated)).multiply(lastBoltzmannWeight);

			ans = ans.add(uniformBound);
		}

		return ans;
	}
	
	
	public double getPercentQPrime() {
		
		BigDecimal ans = BigDecimal.ZERO;
		
		synchronized( lock ) {
			BigInteger remaining = unPruned.subtract(enumerated);
			BigDecimal uniformBound = new BigDecimal(remaining).multiply(lastBoltzmannWeight);
			BigDecimal denom = uniformBound.add(totalPF);
			ans = denom.compareTo(BigDecimal.ZERO) == 0 ? BigDecimal.ONE : totalPF.divide(denom, 4);
		}
		
		return ans.doubleValue();
	}


	protected void nullify() {
		
		if(!confsExhausted) {
			BigInteger remaining = unPruned.subtract(enumerated);
			
			BigDecimal uniformBound = new BigDecimal(remaining).multiply(lastBoltzmannWeight);
			BigDecimal denom = uniformBound.add(totalPF);
			BigDecimal percentQPrime = denom.compareTo(BigDecimal.ZERO) == 0 ? BigDecimal.ONE : totalPF.divide(denom, 4);
			
			System.out.print("q' confTree complete. # enumerated: " + enumerated + ". # remaining: " + remaining + ". ");
			System.out.println("% q': " + percentQPrime);
		}
		confsExhausted = true;
		
		lastBoltzmannWeight = BigDecimal.ZERO;
		pf = null;
		confSearch = null;
	}
}
