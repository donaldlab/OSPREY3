package edu.duke.cs.osprey.kstar.pfunc.impl;

import java.io.Serializable;
import java.math.BigInteger;
import java.util.ArrayList;

import cern.colt.Arrays;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.KSConfQ;
import edu.duke.cs.osprey.kstar.QPrimeCalculator;
import edu.duke.cs.osprey.kstar.RCEnergyContribs;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class PFNew01 extends PFTrad implements Serializable {

	// temp for benchmarking
	protected long startTime;

	protected KSConfQ confsQ = null;
	protected QPrimeCalculator qPrimeCalculator = null;

	public PFNew01() {
		super();
	}

	public PFNew01( int strand, ArrayList<String> sequence, ArrayList<Integer> flexResIndexes, 
			String checkPointPath, String searchProblemName, 
			ConfigFileParser cfp, SearchProblem panSeqSP ) {

		super( strand, sequence, flexResIndexes, checkPointPath, searchProblemName, cfp, panSeqSP );
	}


	public void cleanup() {
		super.cleanup();
	}


	public void start() {

		try {

			setRunState(RunState.STARTED);
			
			if(canUseHotByManualSelection()) 
				createHotsFromCFG();

			// set pstar
			initPStar();

			confsQ = new KSConfQ( this, 1, partialQLB );
			qPrimeCalculator = new QPrimeCalculator( this );
			qPrimeCalculator.setPriority(Thread.MAX_PRIORITY);

			qPrimeCalculator.start();
			confsQ.start();

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}

		startTime = System.currentTimeMillis();
	}


	protected int canSatisfy( int request ) throws Exception {
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

				return granted;
			}
		}

		return request;
	}


	protected void iterate() throws Exception {
		// iterate is only called when eAppx = false
		KSConf conf = null;

		synchronized( confsQ.lock ) {

			if( canSatisfy(1) == 0 )
				return;

			// we are guaranteed that confs can satisfy our request

			// we don't want to hold a lock when we are minimizing, so 
			// we dequeue here and release lock for minimizing
			conf = confsQ.deQueue();

			minimizingConfs = minimizingConfs.add( BigInteger.ONE );

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
		}	
	}


	@Override
	protected void computeSlice() {

		try {

			// this condition only occurs when we are checkpointing
			if( KSAbstract.doCheckPoint) {

				if( !confsQ.isExhausted() && confsQ.getState() == Thread.State.NEW ) {

					// for safety, we can re-start the conformation tree, since i am not
					// entirely sure how cleanly the conformation tree can be serialized and de-serialized
					// confs.restartConfTree();
					confsQ.start();
					synchronized( confsQ.lock ) {
						confsQ.lock.notify();
					}
				}

				if( !qPrimeCalculator.isExhausted() && qPrimeCalculator.getState() == Thread.State.NEW ) {
					qPrimeCalculator.start();
				}
			}

			iterate();

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}

	protected void compute() {

		try {
			// process conformations until e-approximation reached
			while( eAppx == EApproxReached.FALSE ) {

				iterate();

			}
		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	protected void combineResidues(KSConf conf, double pbe, double tpbe, int[] tpce) {

		abort(true);

		System.out.print("% bound error: " + pbe + ". ");
		System.out.print("% bound error from top "+ getHotNumRes() +" RCs: " + tpbe + ". positions: " + Arrays.toString(tpce));		
		System.out.print(". Combining residues: "); for(int i : tpce) System.out.print(getSequence().get(i) + " ");
		System.out.print("... ");

		long start = System.currentTimeMillis();
		sp.mergeResiduePositions(tpce);
		memoizePosInHot(tpce);
		long duration = (System.currentTimeMillis()-start)/1000;

		System.out.println("done in " + duration + "s");
		
		
		//BigDecimal oldPartialQPLB = new BigDecimal(partialQLB.toString());
		partialQLB = reComputePartialQLB(null);
		//if( oldPartialQPLB.compareTo(partialQLB) < 0 )
		//	throw new RuntimeException("ERROR: old partial q' - new partial q': " + oldPartialQPLB.subtract(partialQLB) + " must be >= 0");
		
		
		confsQ = new KSConfQ( this, 1, partialQLB );
		qPrimeCalculator = new QPrimeCalculator( this );
		qPrimeCalculator.setPriority(Thread.MAX_PRIORITY);

		qPrimeCalculator.start();
		confsQ.start();
		
		conf.setEnergyBound(getConfBound(null, conf.getConfArray(), false));
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

		if( isFullyDefined() ) {
			// we do not have a lock when minimizing
			mef = sp.decompMinimizedEnergy(conf.getConfArray());
			energy = mef.getPreCompE();
		}

		else energy = conf.getEnergyBound();

		conf.setEnergy(energy);
		boundError = conf.getEnergyBound() - conf.getEnergy();

		minimizingConfs = minimizingConfs.subtract( BigInteger.ONE );

		updateQStar( conf );

		updateQPrime();

		// negative values of effective esilon are disallowed
		if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0 ) {
			System.out.println("here: " + effectiveEpsilon + " " + qStar + " " + qPrime + " " + pStar);
			eAppx = EApproxReached.NOT_POSSIBLE;
			return;
		}

		long currentTime = System.currentTimeMillis();

		synchronized( confsQ.lock ) {

			if( !PFAbstract.suppressOutput ) {
				if( !printedHeader ) printHeader();
				System.out.println(boundError + "\t" + energy + "\t" + effectiveEpsilon + "\t" + getNumMinimized4Output() + 
						"\t" + getNumUnEnumerated() + "\t" + confsQ.size() + "\t" + ((currentTime-startTime)/1000));
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

		System.out.println("boundError" + "\t" + "minE" + "\t" + "epsilon" + "\t" + "#min" +
				"\t" + "#un-enum" + "\t" + "#buf" + "\t"+ "time(sec)");

		printedHeader = true;
	}


	protected BigInteger getNumUnEnumerated() {
		// assuming locks are in place

		BigInteger numProcessing = getNumMinimized().add(BigInteger.valueOf(confsQ.size())).add(minimizingConfs);

		BigInteger ans = unPrunedConfs.subtract( numProcessing );

		if( ans.compareTo(BigInteger.ZERO) < 0 ) 
			throw new RuntimeException("ERROR: the number of un-enumerated conformations must be >= 0");

		return ans;
	}


	protected void updateQPrime() {

		try {
			while( partialQLB.compareTo(qPrimeCalculator.getTotalQLB()) > 0 ) Thread.sleep(1);
			qPrime = qPrimeCalculator.getQPrime(partialQLB);

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
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

			confsQ.cleanUp(nullify);

			qPrimeCalculator.cleanUp(nullify);

			eAppx = EApproxReached.FALSE;

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	public String getImpl() {
		return "new01";
	}
}
