package edu.duke.cs.osprey.kstar;

import java.math.BigDecimal;
import java.math.BigInteger;

import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;

@SuppressWarnings("serial")
public class QPrimeCalculator extends PStarCalculator {

	private BigInteger unPruned = BigInteger.ZERO;

	public QPrimeCalculator( PFAbstract pf, boolean usePrunedConfs ) {
		super(pf, usePrunedConfs);

		unPruned = pf.getNumUnPruned();
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
