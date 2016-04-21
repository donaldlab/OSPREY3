package edu.duke.cs.osprey.kstar;

import java.io.Serializable;
import java.math.BigDecimal;
import java.math.BigInteger;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;

@SuppressWarnings("serial")
public class QPrimeConfTree extends Thread implements Serializable {

	private PFAbstract pf;
	private BigInteger unPruned;
	private BigDecimal lastBoltzmannWeight = BigDecimal.ZERO;
	private BigInteger enumerated = BigInteger.ZERO;
	private ConfSearch confSearch;
	boolean confsExhausted = false;
	private BigDecimal totalQLB = BigDecimal.ZERO;
	public final String lock = new String("LOCK");

	public QPrimeConfTree( PFAbstract pf, BigInteger unPruned ) {

		this.pf = pf;

		this.unPruned = unPruned.add(BigInteger.ZERO);

		confSearch = pf.getConfTree(false);
	}


	public BigDecimal getQPrime( BigDecimal partialQLB ) {

		BigDecimal ans = BigDecimal.ZERO;

		synchronized( lock ) {
			ans = totalQLB.subtract(partialQLB);

			BigDecimal uniformBound = new BigDecimal(unPruned.subtract(enumerated)).multiply(lastBoltzmannWeight);

			ans = ans.add(uniformBound);
		}

		return ans;
	}


	public BigDecimal getTotalQLB() {
		BigDecimal ans;

		synchronized( lock ) {
			ans = totalQLB;
		}

		return ans;
	}


	public BigInteger getNumEnumerated() {
		BigInteger ans;

		synchronized( lock ) {
			ans = enumerated;
		}

		return ans;
	}


	private void nullify() {
		confsExhausted = true;
		lastBoltzmannWeight = BigDecimal.ZERO;
		pf = null;
		confSearch = null;
		
		System.out.println("qPrimeConfTree complete...");
	}


	public void cleanUp( boolean nullify ) throws InterruptedException {		

		this.join();

		if(nullify) nullify();
	}
	

	public boolean isExhausted() {
		return confsExhausted;
	}


	public void run() {

		if(confsExhausted) { nullify(); return; }

		int conf[];

		while( true ) {

			conf = confSearch.nextConf();

			synchronized( lock ) {

				if( conf == null ) { nullify(); return; }

				double energyBound = pf.getConfBound(confSearch, conf, false);
				if( energyBound == Double.POSITIVE_INFINITY ) { nullify(); return; }

				lastBoltzmannWeight = pf.getBoltzmannWeight(energyBound);
				if( lastBoltzmannWeight.compareTo(BigDecimal.ZERO) == 0 ) { nullify(); return; }

				totalQLB = totalQLB.add( lastBoltzmannWeight );
				enumerated = enumerated.add(BigInteger.ONE);

				// may not be done yet. do not nullify. owner calls cleanup to nullify
				// when we are actually done
				if( pf.getEpsilonStatus() != EApproxReached.FALSE ) { return; }
			}

		}
	}

}
