package edu.duke.cs.osprey.kstar.pfunc.impl;

import java.io.Serializable;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.KSPStarCalculator;
import edu.duke.cs.osprey.kstar.KSQPrimeCalculator;
import edu.duke.cs.osprey.kstar.KSSearchProblem;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class PFParallel0 extends PFParallel1 implements Serializable {
	
	public PFParallel0() {
		super();
	}
	
	
	public PFParallel0( int strand, ArrayList<String> sequence, 
			ArrayList<Integer> absolutePos, 
			String checkPointPath, String reducedSPName, 
			KSConfigFileParser cfp, KSSearchProblem panSP ) {

		super( strand, sequence, absolutePos, checkPointPath, reducedSPName, cfp, panSP );
	}
	
	
	public void start() {

		try {

			setRunState(RunState.STARTED);

			if(canUseHotByManualSelection()) 
				createHotsFromCFG();

			confSearch = getConfTree(false);
			
			// set pstar
			pStarCalculator = null;
			if(prunedConfs.compareTo(BigInteger.ZERO) == 0) pStar = BigDecimal.ZERO;
			else {
				System.out.println("using p* calculator");
				pStarCalculator = new KSPStarCalculator( this );
				pStarCalculator.setPriority(Thread.MAX_PRIORITY);
			}

			qPrimeCalculator = new KSQPrimeCalculator( this );
			qPrimeCalculator.setPriority(Thread.MAX_PRIORITY);

			if(pStarCalculator != null) pStarCalculator.start();
			qPrimeCalculator.start();
			
			if(!isContinuous() && isFullyDefined()) Thread.sleep(initSleepTime);

		} catch (Exception ex) {
			throw new Error("can't compute partition function", ex);
		}

		startTime = System.currentTimeMillis();
	}
	
	
	protected void iterate() {
		try {

			ScoredConf conf;
	
			if( (conf = confSearch.nextConf()) != null ) {
	
				if( processedConfsSet.contains(conf) ) return;
	
				KSConf ksConf = new KSConf(conf.getAssignments(), getConfBound(confSearch, conf.getAssignments()));
	
				accumulate(ksConf);
				
				if( eAppx != EApproxReached.FALSE ) {
					// we leave this function
					qPrimeCalculator.cleanUp(true);
					if(pStarCalculator != null) pStarCalculator.cleanUp(true);
				}
			}
			
			exitIfTimeOut();
			
		} catch (InterruptedException ex) {
			// something interrupted us because it wants us to stop,
			// so throw an exception that no one's supposed to catch
			// and hopefully bubble to the top of the current thread
			throw new Error(ex);
		}
	}
	
	
	protected void accumulate( KSConf conf ) {

		double energy = 0, boundError = 0;
		MultiTermEnergyFunction mef = null;

		if( isContinuous() && isFullyDefined() ) {
			// we do not have a lock when minimizing
			mef = reducedSP.decompMinimizedEnergy(conf.getConfArray());
			energy = mef.getPreCompE();
		}

		else energy = conf.getEnergyBound();

		conf.setEnergy(energy);
		boundError = conf.getEnergyBound() - conf.getEnergy();
		
		updateQStar( conf );

		Et = conf.getEnergyBound();
		updateQPrime();

		// negative values of effective esilon are disallowed
		if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0 ) {

			eAppx = EApproxReached.NOT_POSSIBLE;

			return;
		}

		long currentTime = System.currentTimeMillis();

		if( !PFAbstract.suppressOutput ) {
			if( !printedHeader ) printHeader();

			System.out.println(numberFormat.format(boundError) + "\t" + numberFormat.format(energy) + "\t" 
					+ numberFormat.format(effectiveEpsilon) + "\t" + getNumProcessed() + "\t" 
					+ getNumUnEnumerated() + "\t"+ (currentTime-startTime)/1000);
		}

		eAppx = effectiveEpsilon <= targetEpsilon || maxKSConfsReached() ? EApproxReached.TRUE: EApproxReached.FALSE;

		// hot
		double peb = (conf.getEnergyBound()-conf.getEnergy())/conf.getEnergy();
		if(canUseHotByConfError(peb)) 
			tryHotForConf(conf, mef);

		// for partial sequences when doing KAstar
		if( !isFullyDefined() && eAppx == EApproxReached.TRUE ) adjustQStar();
	}
	
	
	protected void printHeader() {

		System.out.println("error" + "\t" + "energy" + "\t" + "epsilon" + "\t" + "#processed" +
				"\t" + "#un-enum" + "\t" + "time(sec)");

		printedHeader = true;
	}


	public String getImpl() {
		return "Parallel0";
	}


	public void cleanup() {
		super.cleanup();
	}
	
	
	protected BigInteger getNumUnEnumerated() {		
		BigInteger ans = unPrunedConfs.subtract( getNumProcessed() );

		if( ans.compareTo(BigInteger.ZERO) < 0 ) 
			throw new RuntimeException("ERROR: the number of un-enumerated conformations must be >= 0");

		return ans;

	}
	
}
