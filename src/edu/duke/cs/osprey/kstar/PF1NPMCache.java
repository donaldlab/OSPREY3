package edu.duke.cs.osprey.kstar;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class PF1NPMCache extends PFAbstract {

	// temp for benchmarking
	protected long startTime;

	protected KSConfQ confs = null;

	protected PF1NPMCache(ArrayList<String> sequence, ConfigFileParser cfp, 
			SearchProblem sp, PruningControl pc, DEEPerSettings dset, 
			ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones, 
			double EW_I0) {
		super( sequence, cfp, sp, pc, dset, moveableStrands, freeBBZones, EW_I0);
	}


	protected void start() {

		setRunState(RunState.STARTED);

		confs = new KSConfQ( this, (SearchProblem)ObjectIO.deepCopy(sp), 1 );

		// set pstar
		if( confs.getNextConf() ) {

			KSConf conf = confs.peek();

			setPStar( conf.getMinEnergyLowerBound() );
		}

		startTime = System.currentTimeMillis();

		try {

			confs.start();

			if( waitUntilCapacity )
				confs.waitUntilCapacity();

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}

	}


	protected boolean epsilonPossible( int requested ) {
		// eAppx is false when this function is invoked

		// if no confs remaining and not at epsilon, then epsilon is not possible
		if( getNumUnMinimizedConfs().compareTo( BigInteger.valueOf(requested) ) < 0
				|| confs.exhausted() 
				|| confs.getState() == Thread.State.TERMINATED ) { 

			if( eAppx == EApproxReached.FALSE ) {
				eAppx = EApproxReached.NOT_POSSIBLE;
				return false;
			}

		}

		return true;
	}


	protected void iterate() throws Exception {

		KSConf conf = null;

		synchronized( confs.qLock ) {
			// only process when there are confs ready to be processed
			if( confs.size() == 0 ) {
				// eAppx is false when this function is invoked
				epsilonPossible( 1 );	
			}
		}

		if( eAppx == EApproxReached.NOT_POSSIBLE ) {

			confs.cleanUp();

			return;
		}

		synchronized( confs.qLock ) {

			if( !confs.exhausted() && confs.size() == 0 ) confs.qLock.wait();

			// we don't want to hold a lock when we are minimizing, so 
			// we dequeue here and release lock for minimizing
			conf = confs.deQueue();

			// this condition means that confsQ was full (and therefore waiting)
			// before we extracted this conformation, so wake it up
			// it would be wasteful to call notify upon every dequeue operation
			if( confs.getState() == Thread.State.WAITING ) confs.qLock.notify();
		}

		// minimization hapens here
		if( (eAppx = accumulate( conf )) != EApproxReached.FALSE ) {
			// we leave this function
			confs.cleanUp();
		}

		if( eAppx == EApproxReached.NOT_POSSIBLE )
			System.out.println("Cannot reach epsilon");
	}


	@Override
	protected void computeSlice() {

		try {
			/*
			synchronized( confs.qLock ) {
				if( confs.getState() == Thread.State.WAITING ) 
					confs.qLock.notify();
			}
			 */
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


	protected EApproxReached accumulate( KSConf conf ) {
		// we do not have a lock when minimizing	
		conf.setMinEnergy( sp.minimizedEnergy(conf.getConf()) );

		// we need a current snapshot of qDagger, so we lock here
		synchronized( confs.qLock ) {
			// update q*, qDagger, and q' atomically
			updateQStar( conf );

			// update qdagger
			confs.setQDagger( confs.getQDagger().subtract(getBoltzmannWeight(conf.getMinEnergyLowerBound())) );

			Et = confs.size() > 0 ? confs.peekTail().getMinEnergyLowerBound() : conf.getMinEnergyLowerBound();

			// negative values of effective esilon are disallowed
			if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0) return EApproxReached.NOT_POSSIBLE;

			long currentTime = System.currentTimeMillis();

			System.out.println(Et + "\t" + conf.getMinEnergy() + "\t" + effectiveEpsilon + "\t" + 
					getNumMinimizedConfs() + "\t" + getNumUnMinimizedConfs() + "\t "+ ((currentTime-startTime)/1000));
		}

		return effectiveEpsilon > targetEpsilon ? EApproxReached.FALSE : EApproxReached.TRUE;
	}


	public BigDecimal getLowerBound() {

		return getQStar();
	}


	public BigDecimal getUpperBound() {
		// this gives the precise 1-epsilon upperbound for the partition function
		synchronized( confs.qLock ) {

			if( eAppx == EApproxReached.TRUE ) return getQStar();

			updateQPrime();

			BigDecimal uB = (qStar.add(confs.getQDagger())).add(qPrime);

			return uB;
		}
	}


	public BigDecimal getUpperBoundAtEpsilon() {
		// this gives the precise 1-epsilon upperbound for the partition function
		// 1 - (q* + remainder)/(denom) = targetEpsilon

		synchronized( confs.qLock ) {

			if( eAppx == EApproxReached.TRUE ) return getQStar();

			updateQPrime();

			BigDecimal uB;
			BigDecimal denom = ((qStar.add(confs.getQDagger())).add(qPrime)).add(pStar);

			BigDecimal remainder = (BigDecimal.valueOf(1.0 - PFAbstract.targetEpsilon).multiply(denom)).subtract(qStar);

			if( remainder.compareTo( confs.getQDagger().add(qPrime) ) <= 0 )
				uB = qStar.add(remainder);

			else 
				uB = denom.subtract(pStar);

			return uB;
		}
	}


	protected void updateQPrime() {
		
		qPrime = getBoltzmannWeight( Et ).
				multiply( new BigDecimal(getNumUnMinimizedConfs().longValue() 
						- confs.size() - minimizingConfs.longValue() ) );
	}


	protected double computeEffectiveEpsilon() {
		
		updateQPrime();

		BigDecimal divisor = ( (qStar.add(confs.getQDagger())).add(qPrime) ).add(pStar);

		// divisor is 0 iff qstar = 0. this means the energies are too high
		// so epsilon can never be reached
		if( divisor.compareTo(BigDecimal.ZERO) == 0 ) return -1.0;

		BigDecimal maxQStar = divisor.subtract(pStar);

		double minEpsilon = BigDecimal.ONE.subtract( maxQStar.divide(divisor, 4) ).doubleValue();

		if( minEpsilon > targetEpsilon ) return -1.0;

		return BigDecimal.ONE.subtract( qStar.divide(divisor,4) ).doubleValue();
	}


	protected BigInteger getNumUnMinimizedConfs() {
		// assuming locks are in place

		BigInteger numProcessing = (minimizedConfs.add(BigInteger.valueOf(confs.size()))).add(minimizingConfs);

		BigInteger ans = initialUnPrunedConfs.subtract( numProcessing );

		// this final comparison is necessary to maintain the count of remaining
		// confs after we re-start k* due to epsilonnotpossible
		if( ans.compareTo(BigInteger.ZERO) < 0 ) ans = BigInteger.ZERO;

		return ans;
	}


	protected double getStopThreshold() {
		// only valid when num remaining confs > 0
		if( getNumUnMinimizedConfs().compareTo(BigInteger.ZERO) == 0 )
			return Double.POSITIVE_INFINITY;

		return BigDecimal.valueOf(-RT).multiply
				( ( e.log
						( ( qStar.multiply
								( BigDecimal.valueOf(rho) ) ).subtract
								( pStar.subtract
										( confs.getQDagger() ) ) ) ).subtract
						( e.log
								( new BigDecimal(getNumUnMinimizedConfs()) ) ) ).doubleValue();
	}


	public void abort() {
		// k* a* needs this method in order to kill confsq for pruned k* calculations
		// uncertain whether this is the best way to implement this
		try {
			if( getRunState() == RunState.NOTSTARTED ) return;

			// eAppx with true or notpossible have already terminated
			if(eAppx == EApproxReached.FALSE) {

				eAppx = EApproxReached.ABORTED;

				confs.cleanUp();
			}
		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}

}
