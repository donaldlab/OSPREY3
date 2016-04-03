package edu.duke.cs.osprey.kstar.pfunc.impl;

import java.io.Serializable;
import java.math.BigInteger;
import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.KSConfQ;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;


/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class PFNew03 extends PFNew02 implements Serializable {

	ArrayList<MinimizerFiber> slaves = new ArrayList<>();

	protected int fibers;
	protected int threadsPerFiber;
	protected int confsPerThread;
	protected int sleepInterval = 10;

	public PFNew03( int strand, ArrayList<String> sequence, ArrayList<Integer> flexResIndexes, 
			String checkPointPath, String searchProblemName, 
			ConfigFileParser cfp, SearchProblem panSeqSP ) {

		super( strand, sequence, flexResIndexes, checkPointPath, searchProblemName, cfp, panSeqSP );
	}


	public void start() {

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

			ArrayList<SearchProblem> sps = parallelCreateSPs(sp, threadsPerFiber);

			slaves.add( new MinimizerFiber( sps, threadsPerFiber, confsPerThread ) );

			slaves.get(slaves.size()-1).id = i;
		}

		confs = new KSConfQ( this, sp, confsPerThread );

		// set pstar
		setPStar( confs.getNextConfELB() );

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
					if( slave.inputConsumed ) {

						accumulate( slave.buf, true );

						if( eAppx != EApproxReached.FALSE ) {

							slave.inputConsumed = false;

							//System.out.println(slave.id + " initiating termination with " + eAppx);
							break;
						}
					}

					// fill slave conformation buffer
					// lock conformation queue
					synchronized( confs.qLock ) {
						
						int request = slave.buf.size();
						int granted = 0;
						
						if( (granted = canSatisfy(request)) == 0 )
							break;
						
						// reduce the size of buf to match
						while( slave.buf.size() > granted ) {
							slave.buf.remove(slave.buf.size()-1);
						}
						
						for( int i = 0; i < granted; ++i ) {
							slave.buf.set(i, confs.deQueue());
						}
						
						minimizingConfs = minimizingConfs.add( BigInteger.valueOf(slave.buf.size()) );

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

			iterate();

			if( eAppx == EApproxReached.FALSE ) Thread.sleep(sleepInterval);

			else {
				confs.cleanUp(true);

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

			confs.cleanUp(true);

			cleanUpSlaves( eAppx );

			// wait for all slaves to finish
			for( MinimizerFiber slave : slaves ) slave.join();

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	protected void cleanUpSlaves( EApproxReached val ) throws Exception {

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

					accumulate( slave.buf, true );
					
					if( eAppx == EApproxReached.TRUE ) {
					
						val = EApproxReached.TRUE;
					}
				}

				//System.out.println(slave2.id + " buffer notified");
				slave.bufLock.notify();
			}
		}
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

		/**
		 * 
		 * @param confs.size is an integer multiple of this.confs.size
		 * @return
		 */
		public void computeQStar() {

			for( int i = 0; i < buf.size(); ++i ) {

				confs.set(i % sps.size(), buf.get(i));

				if( (i+1) % sps.size() == 0 || (i+1) == buf.size() ) {
					
					// reduce conf size if necessary
					while( confs.size() > buf.size() ) {
						
						confs.remove(confs.size()-1);
						
						indexes.remove(indexes.size()-1);
					}
					
					// sp concurrency reached. update partial partition function
					indexes.parallelStream().forEach( j -> confs.get(j).setMinEnergy(sps.get(j).minimizedEnergy(confs.get(j).getConfArray())) );
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

	
	public static String getImpl() {
		return "new03";
	}
}
