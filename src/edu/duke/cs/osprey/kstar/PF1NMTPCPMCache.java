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
public class PF1NMTPCPMCache extends PF1NPCPMCache {

	ArrayList<MinimizerFiber> slaves = new ArrayList<>();

	protected int fibers;
	protected int threadsPerFiber;
	protected int confsPerThread;
	protected int sleepInterval = 10;

	protected PF1NMTPCPMCache(ArrayList<String> sequence, ConfigFileParser cfp, 
			SearchProblem sp, PruningControl pc, DEEPerSettings dset, 
			ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones, 
			double EW_I0) {

		super( sequence, cfp, sp, pc, dset, moveableStrands, freeBBZones, EW_I0 );
	}


	protected void start() {

		setRunState(RunState.STARTED);

		// make a number of RemoteMinimizers
		fibers = PFAbstract.getNumFibers();
		threadsPerFiber = PFAbstract.getNumThreads();
		confsPerThread = PFAbstract.getConfsThreadBuffer() * threadsPerFiber;

		// #confs sent by master to each slave is a multiple of the number of 
		// sps available at each slave. if #confs in master < the number sent to
		// each slave, the master should just minimize these confs.
		slaves.clear();
		for( int i = 0; i < fibers; ++i ) {

			ArrayList<SearchProblem> sps = new ArrayList<>();
			for( int j = 0; j < threadsPerFiber; ++j ) sps.add((SearchProblem)ObjectIO.deepCopy(sp));

			slaves.add( new MinimizerFiber( sps, threadsPerFiber, confsPerThread ) );

			slaves.get(slaves.size()-1).id = i;
		}

		confs = new KSConfQ( this, sp, confsPerThread );

		// set pstar
		if( confs.getNextConf() ) {

			KSConf conf = confs.peek();

			setPStar( conf.getMinEnergyLowerBound() );
		}

		startTime = System.currentTimeMillis();

		try {

			// start conformation queue
			confs.start();

			if( waitUntilCapacity )
				confs.waitUntilCapacity();

			// start slave threads
			for( MinimizerFiber slave : slaves ) slave.start();

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	protected void iterate() throws Exception {

		for( MinimizerFiber slave : slaves ) {

			if( slave.getState() == Thread.State.WAITING ) {
				// this slave can accept new conformations

				// lock slave
				synchronized( slave.bufLock ) {

					// update global qStar
					if( slave.inputConsumed && (eAppx = accumulate( slave.buf )) != EApproxReached.FALSE ) {

						slave.inputConsumed = false;

						//System.out.println(slave.id + " initiating termination with " + eAppx);
						break;
					}

					// fill slave conformation buffer
					// lock conformation queue
					synchronized( confs.qLock ) {

						// only process when there are enough confs ready to be processed
						if( confs.size() < confsPerThread ) {

							if( !epsilonPossible(confsPerThread) ) {
								//System.out.println(slave.id + " cannot reach epsilon");
								break;
							}

							//System.out.println( slave.id + " waiting to fill buf" );
							if( confs.size() < confsPerThread ) confs.qLock.wait();
						}

						for( int i = 0; i < confsPerThread; ++i ) slave.buf.set(i, confs.deQueue());

						//System.out.println(slave.id + " buf full");

						minimizingConfs = minimizingConfs.add( BigInteger.valueOf(confsPerThread) );

						slave.bufFull = true; slave.inputConsumed = false;

						if( confs.getState() == Thread.State.WAITING ) confs.qLock.notify();
					}

					// notify slave to consume buffer
					slave.bufLock.notify();
				}
			}
		}

	}


	protected void computeSlice() {

		try {
			/*
			synchronized( confs.qLock ) {
				if( confs.getState() == Thread.State.WAITING ) 
					confs.qLock.notify();
			}
			 */
			iterate();

			if( eAppx == EApproxReached.FALSE ) Thread.sleep(sleepInterval);

			else {
				confs.cleanUp();

				cleanUpSlaves( eAppx );

				// wait for all slaves to finish
				for( MinimizerFiber slave : slaves ) slave.join();
			}

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	protected void compute() {

		try {

			while( eAppx == EApproxReached.FALSE ) {

				iterate();

				if( eAppx == EApproxReached.FALSE ) Thread.sleep(sleepInterval);
			}

			confs.cleanUp();

			cleanUpSlaves( eAppx );

			// wait for all slaves to finish
			for( MinimizerFiber slave : slaves ) slave.join();

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	protected void cleanUpSlaves( EApproxReached val ) throws InterruptedException {

		//System.out.println(slave.id + " cleaning up other slaves");

		for( MinimizerFiber slave : slaves ) {

			//System.out.println(slave2.id + " cleaning up");

			slave.setEpsilon(val);

			synchronized( slave.bufLock ) {

				// wait til this thread is no longer runnable
				// there is no notification for this condition, so sleep in a loop
				while( slave.getState() == Thread.State.RUNNABLE ) Thread.sleep(sleepInterval);

				// outstanding remote minimizers can change eappx = false or
				// eappx = notposible to eappx = true
				if( eAppx != EApproxReached.TRUE && slave.inputConsumed ) {

					// if( (eAppx = accumulate2( slave.partialQStar, slave.buf )) == EApproxReached.TRUE ) {
					if( (eAppx = accumulate( slave.buf )) == EApproxReached.TRUE ) {
						val = EApproxReached.TRUE;
					}
				}

				//System.out.println(slave2.id + " buffer notified");
				slave.bufLock.notify();
			}
		}
	}


	protected void updateQStar( BigDecimal partialQ, ArrayList<KSConf> partialQConfs ) {
		qStar = qStar.add( partialQ );

		if(saveTopConfsAsPDB) {
			for( KSConf conf : partialQConfs )
				saveTopConf(conf);
		}

		minimizedConfs = minimizedConfs.add( BigInteger.valueOf(partialQConfs.size()) );
		minimizedConfsDuringInterval = minimizedConfsDuringInterval.add( BigInteger.valueOf(partialQConfs.size()) );
	}


	protected EApproxReached accumulate( BigDecimal partialQ, ArrayList<KSConf> partialQConfs ) {

		// we need a current snapshot of qDagger, so we lock here
		synchronized( confs.qLock ) {
			// update q*, qDagger, minimizingConfs, and q' atomically
			updateQStar( partialQ, partialQConfs );

			minimizingConfs = minimizingConfs.subtract( BigInteger.valueOf(partialQConfs.size()) );

			confs.setQDagger( confs.getQDagger().subtract( computePartialQDagger(partialQConfs) ) );

			Et = confs.peekTail() != null ? confs.peekTail().getMinEnergyLowerBound() 
					: partialQConfs.get(partialQConfs.size()-1).getMinEnergyLowerBound();

			// negative values of effective esilon are disallowed
			if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0) return EApproxReached.NOT_POSSIBLE;

			long currentTime = System.currentTimeMillis();

			System.out.println(partialQConfs.get(partialQConfs.size()-1).getMinEnergy() + "\t" + effectiveEpsilon + "\t" + 
					getNumMinimizedConfs() + "\t" + getNumUnMinimizedConfs() + "\t "+ ((currentTime-startTime)/1000));
		}

		return effectiveEpsilon > targetEpsilon ? EApproxReached.FALSE : EApproxReached.TRUE;
	}
	
	
	@Override
	protected EApproxReached accumulate( ArrayList<KSConf> partialQConfs ) {

		double E = 0;
		
		// we need a current snapshot of qDagger, so we lock here
		synchronized( confs.qLock ) {
			// update q*, qDagger, minimizingConfs, and q' atomically
			confs.setQDagger( confs.getQDagger().subtract( computePartialQDagger(partialQConfs) ) );

			Et = confs.peekTail() != null ? confs.peekTail().getMinEnergyLowerBound() 
					: partialQConfs.get(partialQConfs.size()-1).getMinEnergyLowerBound();

			for( KSConf conf : partialQConfs ) {
				
				minimizingConfs = minimizingConfs.subtract( BigInteger.ONE );
				
				E = conf.getMinEnergy();
				updateQStar( conf );
				
				// negative values of effective epsilon are disallowed
				if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0 ) return EApproxReached.NOT_POSSIBLE;
			
				if( effectiveEpsilon <= targetEpsilon ) break;
			}

			long currentTime = System.currentTimeMillis();

			System.out.println(Et + "\t" + E + "\t" + effectiveEpsilon + "\t" + 
					getNumMinimizedConfs() + "\t" + getNumUnMinimizedConfs() + "\t" + ((currentTime-startTime)/1000));
		}

		return effectiveEpsilon > targetEpsilon ? EApproxReached.FALSE : EApproxReached.TRUE;
	}


	protected class MinimizerFiber extends Thread {

		int id = 0;

		final Object bufLock = new Object();
		ArrayList<SearchProblem> sps = null;
		ArrayList<KSConf> confs = new ArrayList<>();
		ArrayList<Integer> indexes = new ArrayList<>();

		ArrayList<KSConf> buf = new ArrayList<>();
		boolean bufFull = false;
		boolean inputConsumed = false;

		BigDecimal partialQStar = BigDecimal.ZERO;
		EApproxReached eAppx = EApproxReached.FALSE;

		/**
		 * 
		 * @param bufLock governs use of the conformation buffer
		 * @param spsPerMinimizer = number of sps; also the size of the confs list
		 */
		public MinimizerFiber( ArrayList<SearchProblem> sps, 
				int spsPerMinimizer, int confsPerMinimizer ) {

			this.sps = sps;
			sps.trimToSize();

			for( int i = 0; i < spsPerMinimizer; ++i ) { 
				confs.add(null);
				indexes.add(i);
			}
			confs.trimToSize();
			indexes.trimToSize();

			for( int i = 0; i < confsPerMinimizer; ++i ) buf.add(null);
			buf.trimToSize();
		}

		public void setEpsilon( EApproxReached val ) {
			eAppx = val;
			//System.out.println(id + " epsilon set to " + val);
		}

		public BigDecimal getQStar() {
			return partialQStar;
		}

		/**
		 * 
		 * @param confs.size is an integer multiple of this.confs.size
		 * @return
		 */
		public void computeQStar() {

			partialQStar = BigDecimal.ZERO;

			for( int i = 0; i < buf.size(); ++i ) {

				confs.set(i % sps.size(), buf.get(i));

				if( (i+1) % sps.size() == 0 ) {
					// sp concurrency reached. update partial partition function
					indexes.parallelStream().forEach( j -> confs.get(j).setMinEnergy(sps.get(j).minimizedEnergy(confs.get(j).getConf())) );

					for( KSConf conf : confs ) partialQStar =  partialQStar.add( getBoltzmannWeight( conf.getMinEnergy() ) );
				}
			}

			inputConsumed = true;
			bufFull = false;
		}

		public void run() {

			//System.out.println(id + " starting");

			try {

				// until epsilon
				while( eAppx == EApproxReached.FALSE ) {

					//System.out.println(id + " locking");

					synchronized( bufLock ) {

						while( !bufFull && eAppx == EApproxReached.FALSE ) {
							//System.out.println(id + " waiting with eAppx " + eAppx);
							bufLock.wait();
						}

						//System.out.println(id + " resuming with " + eAppx);

						if( eAppx != EApproxReached.FALSE ) {
							break;
						}

						// within this block, we are guaranteed to have access to all
						// confs in the buffer, which has been filled by master
						computeQStar();

						bufLock.notify();
					}
				}

				//System.out.println(id + " exiting");

			} catch(Exception e) {
				System.out.println(e.getMessage());
				e.printStackTrace();
				System.exit(1);
			}
		} 

	}

}
