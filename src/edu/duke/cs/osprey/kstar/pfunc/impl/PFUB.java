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
