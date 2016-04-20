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
	private BigDecimal boltzmannWeight;
	private BigInteger enumerated = BigInteger.ZERO;
	private ConfSearch confSearch;
	private BigDecimal qLBTotal = BigDecimal.ZERO;
	public final String lock = new String("LOCK");

	public QPrimeConfTree( PFAbstract pf, BigInteger unPruned ) {

		this.pf = pf;

		this.unPruned = unPruned.add(BigInteger.ZERO);

		confSearch = pf.getConfTree(false);
	}


	public BigDecimal getQPrime( BigDecimal qLBPartial ) {
		
		BigDecimal ans = BigDecimal.ZERO;

		synchronized( lock ) {

			ans = qLBTotal.subtract(qLBPartial);

			BigDecimal uniformBound = new BigDecimal(unPruned.subtract(enumerated)).multiply(boltzmannWeight);

			ans = ans.add(uniformBound);
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
		pf = null;
		confSearch = null;
	}


	public void cleanup() throws InterruptedException {		
		this.join();
		nullify();
	}


	public void run() {

		int conf[];

		while( true ) {

			synchronized( lock ) {

				conf = confSearch.nextConf();

				if( conf == null || pf.getEpsilonStatus() != EApproxReached.FALSE ) {
					nullify();
					return;
				}

				double energyBound = pf.getConfBound(confSearch, conf, false);

				if( energyBound == Double.POSITIVE_INFINITY ) {
					nullify();
					return;
				}

				boltzmannWeight = pf.getBoltzmannWeight(energyBound);

				if( boltzmannWeight.compareTo(BigDecimal.ZERO) == 0 ) {
					nullify();
					return;
				}

				qLBTotal = qLBTotal.add( boltzmannWeight );

				enumerated = enumerated.add(BigInteger.ONE);
			}
		}
	}

}
