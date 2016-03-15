package edu.duke.cs.osprey.kstar.pfunc.impl;

import java.io.Serializable;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.KSConf;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class PF1NUBNM extends PF1NNoCache implements Serializable {

	protected BigInteger enumeratedConfs = BigInteger.ZERO;

	public PF1NUBNM( ArrayList<String> sequence, String checkPointPath, 
			ConfigFileParser cfp, SearchProblem sp, double EW_I0 ) {

		super( sequence, checkPointPath, cfp, sp, EW_I0 );
	}


	protected void updateQStarUB( KSConf conf ) {

		if(saveTopConfsAsPDB) {
			saveTopConf(conf);
		}

		qStar = qStar.add( getBoltzmannWeight( conf.getMinEnergyLB() ) );
		enumeratedConfs = enumeratedConfs.add(BigInteger.ONE);
	}


	protected BigInteger getNumUnEnumerated() {
		return unPrunedConfs.subtract(enumeratedConfs);
	}


	protected void updateQPrime() {
		qPrime = ( getBoltzmannWeight( Et )).multiply( new BigDecimal(getNumUnEnumerated()) );
	}


	protected double computeDelta() {

		// subsequent boltzmann weights are guaranteed to be lower
		if( qStar.compareTo(BigDecimal.ZERO) == 0 ) return -1.0;
		
		updateQPrime();

		double minDelta = pStar.divide(qStar.add(qPrime), 4).doubleValue();
		if( minDelta > targetEpsilon ) return -1.0;
		
		double delta = (qPrime.add(pStar)).divide(qStar, 4).doubleValue();

		return delta;
	}


	protected void accumulate( KSConf conf ) {

		Et = conf.getMinEnergyLB();

		updateQStarUB( conf );

		if( (effectiveEpsilon = computeDelta()) < 0 ) {
			
			eAppx = EApproxReached.NOT_POSSIBLE;

			return;
		}

		long currentTime = System.currentTimeMillis();

		if( !printedHeader ) printHeader();

		if( enumeratedConfs.longValue() % 2048 == 0 ) {
			System.out.println(Et + "\t" + effectiveEpsilon + "\t" 
					+ enumeratedConfs + "\t" + getNumUnEnumerated() + "\t"+ (currentTime-startTime)/1000);
		}
		
		eAppx = effectiveEpsilon <= targetEpsilon || maxKSConfsReached() ? EApproxReached.TRUE: EApproxReached.FALSE;
		
		if( eAppx == EApproxReached.TRUE ) qStar = qStar.multiply(new BigDecimal(1.0 + effectiveEpsilon));
	}


	protected void printHeader() {

		System.out.println("minELB" + "\t" + "delta" + "\t" + "#enum" +
				"\t" + "#un-enum" + "\t" + "time(sec)");

		printedHeader = true;
	}

}
