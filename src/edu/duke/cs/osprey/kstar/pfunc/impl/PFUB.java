package edu.duke.cs.osprey.kstar.pfunc.impl;

import java.io.Serializable;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;

import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.KSSearchProblem;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class PFUB extends PFTraditional implements Serializable {

	protected BigInteger enumeratedConfs = BigInteger.ZERO;

	public PFUB() {
		super();
	}
	
	public PFUB( int strand, ArrayList<String> sequence, 
			ArrayList<Integer> absolutePos, 
			String checkPointPath, String reducedSPName, 
			KSConfigFileParser cfp, KSSearchProblem panSP ) {

		super( strand, sequence, absolutePos, checkPointPath, reducedSPName, cfp, panSP );
	}


	protected void updateQStarUB( KSConf conf ) {
		qStar = qStar.add( getBoltzmannWeight( conf.getEnergyBound() ) );
		enumeratedConfs = enumeratedConfs.add(BigInteger.ONE);
	}


	protected BigInteger getNumUnEnumerated() {
		return unPrunedConfs.subtract(enumeratedConfs);
	}


	protected void updateQPrime() {
		qPrime = ( getBoltzmannWeight( Et )).multiply( new BigDecimal(getNumUnEnumerated()) );
	}


	protected double computeEffectiveEpsilon() {
		
		BigDecimal dividend = qPrime.add(pStar);
		BigDecimal divisor = qStar.add(dividend);
		
		// energies are too high so epsilon can never be reached
		if( divisor.compareTo(BigDecimal.ZERO) == 0 ) 
			return EPSILON_PHASE_2;
		
		if( qStar.add(qPrime).compareTo(BigDecimal.ZERO) == 0 ) 
			return EPSILON_PHASE_2;

		//double minDelta = pStar.divide(qStar.add(qPrime), 4).doubleValue();
		//if( minDelta > targetEpsilon ) 
		//	return EPSILON_PHASE_2;

		double delta = (qPrime.add(pStar)).divide(qStar, 4).doubleValue();

		return delta;
	}


	protected void accumulate( KSConf conf ) {
		
		updateQStarUB( conf );
		
		Et = conf.getEnergyBound();
		
		updateQPrime();

		if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0 ) {

			eAppx = EApproxReached.NOT_POSSIBLE;

			return;
		}

		long currentTime = System.currentTimeMillis();

		if( !PFAbstract.suppressOutput ) {
			if( !printedHeader ) printHeader();

			if( enumeratedConfs.longValue() % 512 == 0 ) {
				System.out.println(numberFormat.format(Et) + "\t" + numberFormat.format(effectiveEpsilon) + "\t" 
						+ enumeratedConfs + "\t" + getNumUnEnumerated() + "\t"+ (currentTime-startTime)/1000);
			}
		}

		eAppx = effectiveEpsilon <= targetEpsilon || maxKSConfsReached() ? EApproxReached.TRUE: EApproxReached.FALSE;

		if( eAppx == EApproxReached.TRUE ) {
			qStar = qStar.multiply(new BigDecimal(1.0 + effectiveEpsilon));
			// for partial sequences when doing KAstar
			if( !isFullyDefined() ) adjustQStar();
		}
	}


	protected void printHeader() {

		System.out.println("minELB" + "\t" + "delta" + "\t" + "#enum" + "\t" + "#un-enum" + "\t" + "time(sec)");

		printedHeader = true;
	}


	public String getImpl() {
		return "UB";
	}
}
