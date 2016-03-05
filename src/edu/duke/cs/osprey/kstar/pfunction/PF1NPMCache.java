package edu.duke.cs.osprey.kstar.pfunction;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.KSConfQ;
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

		super( sequence, cfp, sp, pc, dset, moveableStrands, freeBBZones, EW_I0 );
	}


	public void start() {

		setRunState(RunState.STARTED);

		confs = new KSConfQ( this, (SearchProblem)ObjectIO.deepCopy(sp), 1 );

		// set pstar
		setPStar( confs.getNextConfELB() );

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


	protected int canSatisfy( int request ) throws Exception {
		int granted = 0;

		// wait for queue to be ready
		while( !confs.canSatisfy(request) ) {

			if( !confs.isExhausted() ) {

				if( confs.getState() == Thread.State.WAITING ) {

					confs.setQCapacity(confs.getQCapacity()+request);

					confs.qLock.notify();
				}

				confs.qLock.wait();
			}

			else {
				granted = confs.size();

				if( granted > 0 )
					return granted;

				eAppx = EApproxReached.NOT_POSSIBLE;

				System.out.println("Cannot reach epsilon");

				confs.cleanUp();

				return granted;
			}
		}

		return request;
	}


	protected void iterate() throws Exception {
		// iterate is only called when eAppx = false
		KSConf conf = null;

		synchronized( confs.qLock ) {

			if( canSatisfy(1) == 0 )
				return;

			// we are guaranteed that confs can satisfy our request

			// we don't want to hold a lock when we are minimizing, so 
			// we dequeue here and release lock for minimizing
			conf = confs.deQueue();

			minimizingConfs = minimizingConfs.add( BigInteger.ONE );

			// this condition means that confsQ was full (and therefore waiting)
			// before we extracted this conformation, so wake it up
			// it would be wasteful to call notify upon every dequeue operation
			if( confs.getState() == Thread.State.WAITING ) confs.qLock.notify();
		}

		// minimization hapens here
		accumulate( conf );

		if( eAppx != EApproxReached.FALSE ) {
			// we leave this function
			confs.cleanUp();
		}	
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


	protected void accumulate( KSConf conf ) {
		
		double E = 0;

		// we do not have a lock when minimizing	
		conf.setMinEnergy( sp.minimizedEnergy(conf.getConf()) );

		// we need a current snapshot of qDagger, so we lock here
		synchronized( confs.qLock ) {

			Et = confs.size() > 0 ? confs.peekTail().getMinEnergyLB() : conf.getMinEnergyLB();

			minimizingConfs = minimizingConfs.subtract( BigInteger.ONE );

			// update q*, qDagger, and q' atomically
			E = conf.getMinEnergy();
			updateQStar( conf );

			// update qdagger
			confs.setQDagger( confs.getQDagger().subtract(getBoltzmannWeight(conf.getMinEnergyLB())) );
			if(PFAbstract.useRigEnergy) {
				confs.setQDot( confs.getQDot().subtract( getBoltzmannWeight(conf.getMinEnergyUB()) ) );
			}

			updateQPrime();

			// negative values of effective esilon are disallowed
			if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0 ) {

				eAppx = EApproxReached.NOT_POSSIBLE;

				return;
			}

			long currentTime = System.currentTimeMillis();

			if( !printedHeader ) printHeader();
			
			System.out.println(E + "\t" + effectiveEpsilon + "\t" + 
					getNumMinimizedConfs() + "\t" + getNumUnMinimizedConfs() + "\t" + confs.size() + "\t" + ((currentTime-startTime)/1000));

			eAppx = effectiveEpsilon > targetEpsilon ? EApproxReached.FALSE : EApproxReached.TRUE;
		}

	}


	protected void printHeader() {

		System.out.println("minE" + "\t\t\t" + "epsilon" + "\t\t" + "#min" +
				"\t" + "#un-min" + "\t" + "#buf" + "\t"+ "time(sec)");

		printedHeader = true;
	}


	protected void updateQPrime() {

		qPrime = getBoltzmannWeight( Et ).
				multiply( new BigDecimal(getNumUnMinimizedConfs().longValue() 
						- confs.size() - minimizingConfs.longValue() ) );
	}


	protected double computeEffectiveEpsilon() {

		BigDecimal qPrimePStar = qPrime.add(pStar);

		if(qPrimePStar.compareTo(confs.getCapacityThresh().multiply(confs.getQDagger())) < 0)
			confs.setQCapacity(confs.size());

		else
			confs.restoreQCapacity();

		BigDecimal divisor = ( (qStar.add(confs.getQDagger())).add(qPrimePStar) );

		// divisor is 0 iff qstar = 0. this means the energies are too high
		// so epsilon can never be reached
		if( divisor.compareTo(BigDecimal.ZERO) == 0 ) return -1.0;

		BigDecimal maxQStar = divisor.subtract(pStar);

		double minEpsilon = BigDecimal.ONE.subtract( maxQStar.divide(divisor, 4) ).doubleValue();

		if( minEpsilon > targetEpsilon ) return -1.0;

		BigDecimal dividend = qStar;
		if(PFAbstract.useRigEnergy) {
			dividend = dividend.add( confs.getQDot() );
		}

		return BigDecimal.ONE.subtract( dividend.divide(divisor,4) ).doubleValue();
	}


	protected BigInteger getNumUnMinimizedConfs() {
		// assuming locks are in place

		BigInteger numProcessing = (getNumMinimizedConfs().
				add(BigInteger.valueOf(confs.size()))).add(minimizingConfs);

		BigInteger ans = initialUnPrunedConfs.subtract( numProcessing );

		// this final comparison is necessary to maintain the count of remaining
		// confs after we re-start k* due to epsilonnotpossible
		if( ans.compareTo(BigInteger.ZERO) < 0 ) ans = BigInteger.ZERO;

		return ans;
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
