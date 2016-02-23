package edu.duke.cs.osprey.kstar;

import java.math.BigDecimal;
import java.util.ArrayList;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.kstar.KSConfQ;
import edu.duke.cs.osprey.pruning.PruningControl;


/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class PF1NAStar extends PFAbstract {

	protected long startTime;

	protected KSConfQ confs = null;

	protected PF1NAStar(ArrayList<String> sequence, ConfigFileParser cfp, 
			SearchProblem sp, PruningControl pc, DEEPerSettings dset, 
			ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones, 
			double EW_I0) {
		
		super( sequence, cfp, sp, pc, dset, moveableStrands, freeBBZones, EW_I0);
	}


	protected void start() {

		setRunState(RunState.STARTED);
		
		confs = new KSConfQ( this, sp, 1 );

		// set pstar
		if( confs.getNextConf() ) {

			KSConf conf = confs.peek();

			setPStar( conf.getMinEnergyLowerBound() );
		}
		
		startTime = System.currentTimeMillis();
	}


	protected void computeSlice() {
		try {

			iterate();

		} catch (Exception e) {

			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);

		}
	}


	protected void compute() {

		computeSlice();
	}


	public BigDecimal getQDagger() {

		if(confs == null) 
			return null; 

		synchronized( confs.qLock ) {
			return confs.getQDagger();
		}
	}


	protected void iterate() throws Exception {

		boolean nextConf = false;
		while( (nextConf = confs.getNextConf()) == true && (eAppx = accumulate()) == EApproxReached.FALSE );

		if( !nextConf && eAppx == EApproxReached.FALSE )
			eAppx = EApproxReached.NOT_POSSIBLE;
	}


	protected void updateQPrime() {

		Et = confs.size() > 0 ? confs.get(confs.size()-1).getMinEnergyLowerBound() : 0.0;

		qPrime = getBoltzmannWeight( Et ).
				multiply( new BigDecimal(getNumUnMinimizedConfs().longValue() - confs.size() ) );
	}


	public BigDecimal getUpperBound() {

		return confs.getQDagger();
	}


	public BigDecimal getUpperBound2() {

		synchronized( confs.qLock ) {

			updateQPrime();

			BigDecimal uB = (confs.getQDagger()).add(qPrime);

			return uB;
		}
	}


	public BigDecimal getQStar() {
		return getUpperBound();
	}
	
	
	protected double computeEffectiveEpsilon() {

		updateQPrime();

		BigDecimal divisor = ( qPrime.add(confs.getQDagger()) ).add(pStar);

		// divisor is 0 iff energies are too high, so epsilon can never be reached
		if( divisor.compareTo(BigDecimal.ZERO) == 0 ) return -1.0;

		BigDecimal maxQStar = divisor.subtract(pStar);

		BigDecimal minEpsilon = BigDecimal.ONE.subtract( maxQStar.divide(divisor, 4) );

		if( minEpsilon.compareTo(BigDecimal.valueOf(targetEpsilon)) > 0 ) 
			return -1.0;

		BigDecimal epsilon = BigDecimal.ONE.subtract( confs.getQDagger().divide(divisor, 4) );
		
		return epsilon.doubleValue();
	}


	protected EApproxReached accumulate() {
		// negative values of effective esilon are disallowed
		if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0) return EApproxReached.NOT_POSSIBLE;

		long currentTime = System.currentTimeMillis();

		if( effectiveEpsilon > targetEpsilon )
			return EApproxReached.FALSE;
		
		else {
			System.out.println(effectiveEpsilon + "\t" + 
					getNumMinimizedConfs() + "\t" + getNumUnMinimizedConfs() + "\t "+ ((currentTime-startTime)/1000));
			
			return EApproxReached.TRUE;
		}
	}
	
}
