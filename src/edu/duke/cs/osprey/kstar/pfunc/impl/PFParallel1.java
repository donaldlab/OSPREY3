package edu.duke.cs.osprey.kstar.pfunc.impl;

import java.io.Serializable;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;

import cern.colt.Arrays;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.KSConfQ;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.KSSearchProblem;
import edu.duke.cs.osprey.kstar.KSPStarCalculator;
import edu.duke.cs.osprey.kstar.KSQPrimeCalculator;
import edu.duke.cs.osprey.kstar.RCEnergyContribs;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class PFParallel1 extends PFTraditional implements Serializable {

	protected KSConfQ confsQ = null;
	protected KSQPrimeCalculator qPrimeCalculator = null;
	protected KSPStarCalculator pStarCalculator = null;
	protected long initSleepTime = 3000;

	public PFParallel1() {
		super();
	}

	public PFParallel1( int strand, ArrayList<String> sequence, 
			ArrayList<Integer> absolutePos, 
			String checkPointPath, String reducedSPName, 
			KSConfigFileParser cfp, KSSearchProblem panSP ) {

		super( strand, sequence, absolutePos, checkPointPath, reducedSPName, cfp, panSP );
	}


	public void cleanup() {
		super.cleanup();
		confsQ = null;
		qPrimeCalculator = null;
		pStarCalculator = null;
	}


	public void start() {

		try {

			setRunState(RunState.STARTED);

			if(canUseHotByManualSelection()) 
				createHotsFromCFG();

			// set pstar
			pStarCalculator = null;
			if(prunedConfs.compareTo(BigInteger.ZERO) == 0) pStar = BigDecimal.ZERO;
			else {
				System.out.println("using p* calculator");
				pStarCalculator = new KSPStarCalculator( this );
				pStarCalculator.setPriority(Thread.MAX_PRIORITY);
			}

			confsQ = new KSConfQ( this, 1, partialQLB );
			qPrimeCalculator = new KSQPrimeCalculator( this );
			qPrimeCalculator.setPriority(Thread.MAX_PRIORITY);

			if(pStarCalculator != null) pStarCalculator.start();
			qPrimeCalculator.start();
			confsQ.start();
			
			if(!isContinuous() && isFullyDefined()) Thread.sleep(initSleepTime);

		} catch(Exception ex) {
			throw new Error("can't compute partition function", ex);
		}

		startTime = System.currentTimeMillis();
	}


	protected int canSatisfy( int request ) throws InterruptedException {
		int granted = 0;

		// wait for queue to be ready
		while( !confsQ.canSatisfy(request) ) {

			if( !confsQ.isExhausted() ) {

				if( confsQ.getState() == Thread.State.WAITING || confsQ.getState() == Thread.State.BLOCKED ) {

					if( confsQ.getQCapacity() < request )
						confsQ.setQCapacity(request);

					confsQ.lock.notify();
				}

				confsQ.lock.wait();
			}

			else {

				granted = confsQ.size();

				if( granted > 0 )
					return granted;

				eAppx = EApproxReached.NOT_POSSIBLE;

				// System.out.println("Cannot reach epsilon");

				confsQ.cleanUp(true);
				qPrimeCalculator.cleanUp(true);
				if(pStarCalculator != null) pStarCalculator.cleanUp(true);

				return granted;
			}
		}

		return request;
	}


	protected void iterate() {
		try {
			
			// iterate is only called when eAppx = false
			KSConf conf = null;

			synchronized( confsQ.lock ) {

				if( canSatisfy(1) == 0 )
					return;

				// we are guaranteed that confs can satisfy our request

				// we don't want to hold a lock when we are minimizing, so 
				// we dequeue here and release lock for minimizing
				conf = confsQ.deQueue();

				processingConfs = processingConfs.add( BigInteger.ONE );

				// this condition means that confsQ was full (and therefore waiting)
				// before we extracted this conformation, so wake it up
				// it would be wasteful to call notify upon every dequeue operation
				if( confsQ.getState() == Thread.State.WAITING ) confsQ.lock.notify();
			}

			// minimization hapens here
			accumulate( conf );

			if( eAppx != EApproxReached.FALSE ) {
				// we leave this function
				confsQ.cleanUp(true);
				qPrimeCalculator.cleanUp(true);
				if(pStarCalculator != null) pStarCalculator.cleanUp(true);
			}
			
			exitIfTimeOut();
		
		} catch (InterruptedException ex) {
			// something interrupted us because it wants us to stop,
			// so throw an exception that no one's supposed to catch
			// and hopefully bubble to the top of the current thread
			throw new Error(ex);
		}
	}


	@Override
	protected void computeSlice() {

		try {

			// this condition only occurs when we are checkpointing
			if( KSAbstract.doCheckPoint) {

				if( confsQ != null && !confsQ.isExhausted() && confsQ.getState() == Thread.State.NEW ) {

					// for safety, we can re-start the conformation tree, since i am not
					// entirely sure how cleanly the conformation tree can be serialized and de-serialized
					// confs.restartConfTree();
					confsQ.start();
					synchronized( confsQ.lock ) {
						confsQ.lock.notify();
					}
				}

				if( qPrimeCalculator!= null && !qPrimeCalculator.isExhausted() && qPrimeCalculator.getState() == Thread.State.NEW ) {
					qPrimeCalculator.start();
				}

				if( pStarCalculator != null && !pStarCalculator.isExhausted() && pStarCalculator.getState() == Thread.State.NEW ) {
					pStarCalculator.start();
				}
			}

			iterate();

		} catch(Exception ex) {
			throw new Error("can't compute partition function slice", ex);
		}
	}

	protected void compute() {

		try {
			// process conformations until e-approximation reached
			while( eAppx == EApproxReached.FALSE ) {

				iterate();

			}
		} catch(Exception ex) {
			throw new Error("can't compute partition function", ex);
		}
	}


	protected void combineResidues(KSConf conf, double pbe, double tpbe, int[] tpce) {

		abort(true);

		System.out.print("% bound error: " + pbe + ". ");
		System.out.print("% bound error from top "+ getHotNumRes() +" RCs: " + tpbe + ". positions: " + Arrays.toString(tpce));		
		System.out.print(". Combining residues: "); for(int i : tpce) System.out.print(getSequence().get(i) + " ");
		System.out.print("... ");

		long start = System.currentTimeMillis();
		reducedSP.mergeResiduePositions(tpce);
		memoizePosInHot(tpce);
		long duration = (System.currentTimeMillis()-start)/1000;

		System.out.println("done in " + duration + "s");


		//BigDecimal oldPartialQPLB = new BigDecimal(partialQLB.toString());
		partialQLB = reComputePartialQLB(null);
		//if( oldPartialQPLB.compareTo(partialQLB) < 0 )
		//	throw new RuntimeException("ERROR: old partial q' - new partial q': " + oldPartialQPLB.subtract(partialQLB) + " must be >= 0");


		confsQ = new KSConfQ( this, 1, partialQLB );
		qPrimeCalculator = new KSQPrimeCalculator( this );
		qPrimeCalculator.setPriority(Thread.MAX_PRIORITY);

		qPrimeCalculator.start();
		confsQ.start();

		conf.setEnergyBound(getConfBound(null, conf.getConfArray()));
	}


	protected void tryHotForConf(KSConf conf, MultiTermEnergyFunction mef) {

		double pbe = 0, tpbe = 0; int[] tpce = null;

		RCEnergyContribs rce = new RCEnergyContribs(this, mef, conf.getConfArray());
		pbe = rce.getPercentBoundError();

		if(pbe >= getHotBoundPct()) {
			tpbe = rce.getPercentErrorForTopPos(getHotNumRes());
			if(tpbe >= getHotTopRotsPct()) {
				tpce = rce.getTopPosCausingError(getHotNumRes());
				combineResidues(conf, pbe, tpbe, tpce);
			}
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

		processingConfs = processingConfs.subtract( BigInteger.ONE );

		updateQStar( conf );

		updateQPrime();

		updatePStar();

		// negative values of effective esilon are disallowed
		if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0 ) {
			eAppx = EApproxReached.NOT_POSSIBLE;
			return;
		}

		long currentTime = System.currentTimeMillis();

		synchronized( confsQ.lock ) {

			if( !PFAbstract.suppressOutput ) {
				if( !printedHeader ) printHeader();

				System.out.println(numberFormat.format(boundError) + "\t" + numberFormat.format(energy) + "\t" 
						+ numberFormat.format(effectiveEpsilon) + "\t" + getNumProcessed() + "\t" 
						+ getNumUnEnumerated() + "\t" + confsQ.size() + "\t" + ((currentTime-startTime)/1000));
			}
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
				"\t" + "#un-enum" + "\t" + "#buf" + "\t"+ "time(sec)");

		printedHeader = true;
	}


	protected BigInteger getNumUnEnumerated() {
		// assuming locks are in place
		
		BigInteger numProcessing = getNumProcessed().add(BigInteger.valueOf(confsQ.size())).add(processingConfs);

		BigInteger ans = unPrunedConfs.subtract( numProcessing );

		if( ans.compareTo(BigInteger.ZERO) < 0 ) 
			throw new RuntimeException("ERROR: the number of un-enumerated conformations must be >= 0");

		return ans;
	}


	protected void updateQPrime() {

		try {

			if(qPrimeCalculator == null) return;

			while( partialQLB.compareTo(qPrimeCalculator.getTotalPF()) > 0 ) Thread.sleep(1000);

			qPrime = qPrimeCalculator.getQPrime(partialQLB);

		} catch(Exception ex) {
			throw new Error("can't compute partition function q'", ex);
		}
	}


	protected void updatePStar() {

		try {

			if(pStarCalculator == null) return;

			pStar = pStarCalculator.getPStar();

		} catch(Exception ex) {
			throw new Error("can't compute partition function p*", ex);
		}
	}


	public void abort(boolean nullify) {
		// k* a* needs this method in order to kill confsq for pruned k* calculations
		// uncertain whether this is the best way to implement this
		try {
			if( getRunState() == RunState.NOTSTARTED ) return;

			// other values of eAppx have already terminated
			if(eAppx != EApproxReached.FALSE) return;

			eAppx = EApproxReached.NOT_POSSIBLE;

			if(confsQ != null) confsQ.cleanUp(nullify);

			if(qPrimeCalculator != null) qPrimeCalculator.cleanUp(nullify);

			if(pStarCalculator != null) pStarCalculator.cleanUp(nullify);

			eAppx = EApproxReached.FALSE;

		} catch(Exception ex) {
			throw new Error("can't abort partition function computation", ex);
		}
	}


	public String getImpl() {
		return "Parallel1";
	}
}
