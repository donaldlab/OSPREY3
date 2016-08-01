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
