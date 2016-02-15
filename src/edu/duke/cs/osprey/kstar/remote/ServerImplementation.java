package edu.duke.cs.osprey.kstar.remote;

import java.rmi.RemoteException;
import java.rmi.server.UnicastRemoteObject;
import java.util.ArrayList;
import java.util.Collections;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.minimization.MinimizerFactory;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
import edu.duke.cs.osprey.tools.ObjectIO;


@SuppressWarnings("serial")
public class ServerImplementation extends UnicastRemoteObject implements ServerInterface {

	boolean envVarsSet = false;
	
	ArrayList<RemoteMinimizer> slaves = new ArrayList<>();
	int fibers = 0;
	int confsPerThread = 0;
	int threadsPerFiber = 0;

	// increments from 0 to numMinimizers. finished when numfinished = numminimizers
	int numFinished = 0;
	boolean processed = false;
	final Object processedLock = new Object();

	boolean ready = true;
	boolean epsilon = false;

	int processedConfsSize = 0;
	ArrayList<KSConf> processedConfs = new ArrayList<>();

	protected ServerImplementation() throws RemoteException {
		super();
		// TODO Auto-generated constructor stub
	}

	
	public boolean envVarsSet() throws RemoteException {
		return envVarsSet;
	}
	

	public void setEnvVars( int fibers, int threadsPerFiber, 
			int confsPerThread,
			String eMinMethod,
			EnergyFunctionGenerator curEFcnGenerator,
			GenericResidueTemplateLibrary resTemplates,
			boolean assignTemplatesToStruct,
			boolean deleteNonTemplateResidues,
			boolean useMPI,
			double DUNBRACK_PROBABILTY_CUTOFF,
			String dataDir ) throws RemoteException {

		// clear previous stuff
		softTerminate();
		slaves.clear();
		
		EnvironmentVars.curEFcnGenerator = curEFcnGenerator;
		EnvironmentVars.resTemplates = resTemplates;
		EnvironmentVars.assignTemplatesToStruct = assignTemplatesToStruct;
		EnvironmentVars.deleteNonTemplateResidues = deleteNonTemplateResidues;
		EnvironmentVars.useMPI = useMPI;
		EnvironmentVars.DUNBRACK_PROBABILTY_CUTOFF = DUNBRACK_PROBABILTY_CUTOFF;
		EnvironmentVars.setDataDir(dataDir);

		this.fibers = fibers;
		this.confsPerThread = confsPerThread;
		this.threadsPerFiber = threadsPerFiber;

		MinimizerFactory.setImplementation( eMinMethod );

		// init common forkjoinpool to #hardware threads-1
		int threads = Runtime.getRuntime().availableProcessors()-1;
		ThreadParallelism.setNumThreads( threads );
		
		envVarsSet = true;
	}

	@Override
	public void initFibers( SearchProblem sp ) throws RemoteException {

		// clear previous stuff
		softTerminate();
		slaves.clear();

		// initialize processed confs arraylist
		processedConfs.clear();
		for( int i = 0; i < fibers * confsPerThread; ++i ) processedConfs.add(null);
		processedConfs.trimToSize();

		for( int i = 0; i < fibers; ++i ) {

			ArrayList<SearchProblem> sps = new ArrayList<>();
			for( int j = 0; j < threadsPerFiber; ++j ) sps.add((SearchProblem)ObjectIO.deepCopy(sp));

			RemoteMinimizer slave = new RemoteMinimizer( this, sps, threadsPerFiber, confsPerThread );
			slaves.add( slave );

			slave.start();
		}
	}


	@Override
	public void setUnprocessedConfs( ArrayList<KSConf> unProcessedConfs ) throws RemoteException {

		ready = false;

		// check that size of unProcessedConfs is numMinimizers * confsPerMinimizer
		if( unProcessedConfs.size() != fibers * confsPerThread )
			throw new RuntimeException("ERROR: acceptUnprocessedConfs did not receive "
					+ "numMinimizers * confsPerMinimizer = " + fibers * confsPerThread
					+ " conformations.");

		// assign confsPerMinimizer confs to each remote minimizer
		int i = 0;
		for( RemoteMinimizer slave : slaves ) {

			for( int j = 0; j < confsPerThread; ++j ) {
				slave.buf.set(j, unProcessedConfs.get(i++) );
			}

			synchronized( slave.waitLock ) {
				slave.waitLock.notify();
			}

		}
	}


	@SuppressWarnings("unchecked")
	@Override
	public ArrayList<KSConf> getProcessedConfs() throws RemoteException {

		// this is needed for correct q''
		Collections.sort( processedConfs );

		// clean up state
		processedConfsSize = 0;
		numFinished = 0;

		processed = false;
		ready = true;

		return processedConfs;

	}


	@Override
	public boolean isProcessed() throws RemoteException {

		synchronized( processedLock ) {

			if( numFinished == fibers ) {
				processed = true;
			}

			return processed;
		}

	}


	@Override
	public boolean isReady() throws RemoteException {

		synchronized( processedLock ) {
			return ready;
		}
	}


	@Override
	public void setEpsilon( boolean epsilon ) throws RemoteException {
		this.epsilon = epsilon;
	}


	@Override
	public void softTerminate() throws RemoteException {
		try {

			setEpsilon(true);

			// wait for remote minimizer to finish
			for( RemoteMinimizer slave : slaves ) {

				if( slave.getState() == Thread.State.WAITING ) {

					synchronized( slave.waitLock ) {
						slave.waitLock.notify();
					}
				}

				slave.join();
			}

			processedConfsSize = 0;
			numFinished = 0;

			processed = false;
			ready = true;
			setEpsilon(false);

		} catch (Exception e) {
			System.out.println( e.getMessage() );
			e.printStackTrace();
			System.exit(1);
		}	
	}


	@Override
	public void hardTerminate() throws RemoteException {

		softTerminate();

		System.out.println("Server:\t" + Constants.RMI_ID + " is exiting.");

		System.exit(0);

	}


	protected class RemoteMinimizer extends Thread {

		ServerImplementation parent = null;

		ArrayList<SearchProblem> sps = null;

		ArrayList<KSConf> buf = new ArrayList<>();
		ArrayList<KSConf> confs = new ArrayList<>();
		ArrayList<Integer> indexes = new ArrayList<>();

		final Object waitLock = new Object();

		public RemoteMinimizer( ServerImplementation parent, ArrayList<SearchProblem> sps, 
				int threadsPerFiber, int confsPerThread ) {

			this.parent = parent;

			this.sps = sps;
			sps.trimToSize();

			for( int i = 0; i < threadsPerFiber; ++i ) { 
				confs.add(null);
				indexes.add(i);
			}
			confs.trimToSize();
			indexes.trimToSize();

			for( int i = 0; i < confsPerThread; ++i ) buf.add(null);
			buf.trimToSize();
		}


		public void computeMinimizedEnergies() {

			for( int i = 0; i < buf.size(); ++i ) {

				confs.set(i % sps.size(), buf.get(i));

				if( (i+1) % sps.size() == 0 ) {
					// store minimized energies
					indexes.parallelStream().forEach( j -> confs.get(j).setMinEnergy(sps.get(j).minimizedEnergy(confs.get(j).getConf())) );
				}
			}

			// update number of finished minimizers
			// assemble all processed conformations into the output structure
			synchronized( parent.processedLock ) {

				for( KSConf conf : buf )
					parent.processedConfs.set(parent.processedConfsSize++, conf);

				parent.numFinished++;

			}

		}


		public void run() {
			try {

				while( !parent.epsilon ) {

					// we are either waiting or computing
					synchronized( waitLock ) {

						if( !parent.epsilon ) waitLock.wait();

						if( parent.epsilon ) break;
					}

					computeMinimizedEnergies();
				}

			} catch( Exception e ) {
				System.out.println(e.getMessage());
				e.printStackTrace();
				System.exit(1);
			}
		}

	}

}
