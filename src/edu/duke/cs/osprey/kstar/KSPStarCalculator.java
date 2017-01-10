package edu.duke.cs.osprey.kstar;

import java.io.Serializable;
import java.math.BigDecimal;
import java.math.BigInteger;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;

@SuppressWarnings("serial")
public class KSPStarCalculator extends Thread implements Serializable {

	protected PFAbstract pf; 
	private BigInteger pruned = BigInteger.ZERO;
	protected BigDecimal lastBoltzmannWeight = BigDecimal.ZERO;
	protected double lastEnergyBound = Double.NEGATIVE_INFINITY;
	protected BigInteger enumerated = BigInteger.ZERO;
	protected ConfSearch confSearch = null;
	protected boolean confsExhausted = false;
	protected BigDecimal totalPF = BigDecimal.ZERO;
	protected final String lock = new String("LOCK");

	public KSPStarCalculator( PFAbstract pf ) {

		this.pf = pf;

		pruned = pf.getNumPruned();
		
		confSearch = getConfSearch();
	}
	
	
	private ConfSearch getConfSearch() {
		return pf.getConfTree(true);
	}
	
	
	public BigDecimal getPStar() {

		BigDecimal ans = BigDecimal.ZERO;

		synchronized( lock ) {
			ans = totalPF;
			
			BigDecimal uniformBound = new BigDecimal(pruned.subtract(enumerated)).multiply(lastBoltzmannWeight);

			ans = ans.add(uniformBound);
		}

		return ans;
	}
	
	
	public BigDecimal getTotalPF() {
		BigDecimal ans;

		synchronized( lock ) {
			ans = totalPF;
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
	
	
	protected void nullify() {
		
		if(!confsExhausted) {
			BigInteger remaining = pruned.subtract(enumerated);
			
			BigDecimal uniformBound = new BigDecimal(remaining).multiply(lastBoltzmannWeight);
			BigDecimal denom = uniformBound.add(totalPF);
			BigDecimal percentPStar = denom.compareTo(BigDecimal.ZERO) == 0 ? BigDecimal.ONE : totalPF.divide(denom, 4);
			
			System.out.print("p* confTree complete. # enumerated: " + enumerated + ". # remaining: " + remaining + ". ");
			System.out.println("% p*: " + percentPStar);
		}
		confsExhausted = true;
		
		lastBoltzmannWeight = BigDecimal.ZERO;
		pf = null;
		confSearch = null;
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

		ScoredConf conf;

		while( true ) {

			conf = confSearch.nextConf();

			synchronized( lock ) {

				if( conf == null ) { nullify(); return; }

				lastEnergyBound = pf.getConfBound(confSearch, conf.getAssignments());
				if( lastEnergyBound == Double.POSITIVE_INFINITY ) { 
					lastBoltzmannWeight = BigDecimal.ZERO;
					totalPF = totalPF.add( BigDecimal.ZERO ); enumerated = enumerated.add(BigInteger.ONE);
					nullify(); return; 
				}

				lastBoltzmannWeight = pf.getBoltzmannWeight(lastEnergyBound);
				if( lastBoltzmannWeight.compareTo(BigDecimal.ZERO) == 0 ) { 
					totalPF = totalPF.add( BigDecimal.ZERO ); enumerated = enumerated.add(BigInteger.ONE);
					nullify(); return; 
				}

				totalPF = totalPF.add( lastBoltzmannWeight );
				enumerated = enumerated.add(BigInteger.ONE);

				// may not be done yet. do not nullify. owner calls cleanup to nullify when we are actually done
				// if( pf.getEpsilonStatus() != EApproxReached.FALSE ) { return; }
				if( pf.getEpsilonStatus() != EApproxReached.FALSE ) { nullify(); return; }
			}

		}
	}
}
