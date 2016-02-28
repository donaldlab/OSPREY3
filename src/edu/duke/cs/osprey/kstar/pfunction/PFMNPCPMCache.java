package edu.duke.cs.osprey.kstar.pfunction;

import java.math.BigInteger;
import java.rmi.RemoteException;
import java.rmi.registry.LocateRegistry;
import java.rmi.registry.Registry;
import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.kstar.KSConfQ;
import edu.duke.cs.osprey.kstar.remote.Constants;
import edu.duke.cs.osprey.kstar.remote.ServerInterface;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.pruning.PruningControl;

public class PFMNPCPMCache extends PF1NMTPCPMCache {

	ArrayList<ServerInterface> serverInterfaces = new ArrayList<>();
	ArrayList<KSConf> unProcessedConfs = new ArrayList<>();


	public PFMNPCPMCache(ArrayList<String> sequence, ConfigFileParser cfp, 
			SearchProblem sp, PruningControl pc, DEEPerSettings dset, 
			ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones, 
			double EW_I0) {

		super( sequence, cfp, sp, pc, dset, moveableStrands, freeBBZones, EW_I0);
	}

	public void start() {

		try {

			setRunState(RunState.STARTED);

			fibers = PFAbstract.getNumFibers();
			threadsPerFiber = PFAbstract.getNumThreads();
			confsPerThread = PFAbstract.getConfsThreadBuffer() * threadsPerFiber;

			// establish connection to minimizaiton servers
			serverInterfaces.clear();
			for( String server : serverList ) {

				Registry registry = LocateRegistry.getRegistry( server, Constants.RMI_PORT );

				ServerInterface serverInterface = (ServerInterface)registry.lookup( Constants.RMI_ID );
				serverInterfaces.add( serverInterface );

				if( !serverInterface.envVarsSet() ) {

					System.out.println("Setting environment variables for server: " + server);

					serverInterface.setEnvVars(fibers, threadsPerFiber, confsPerThread,
							PFAbstract.eMinMethod, 
							EnvironmentVars.curEFcnGenerator,
							EnvironmentVars.resTemplates,
							EnvironmentVars.assignTemplatesToStruct,
							EnvironmentVars.deleteNonTemplateResidues,
							EnvironmentVars.useMPI,
							EnvironmentVars.DUNBRACK_PROBABILTY_CUTOFF,
							EnvironmentVars.getDataDir());
				}
			}

			for( ServerInterface serverInterface : serverInterfaces ) {
				serverInterface.initFibers(sp);
			}

			unProcessedConfs.clear();
			for( int i = 0; i < fibers * confsPerThread; ++i ) unProcessedConfs.add(null);
			unProcessedConfs.trimToSize();

			confs = new KSConfQ( this, sp, unProcessedConfs.size() );

			// set pstar
			if( confs.getNextConf() ) {

				KSConf conf = confs.peek();

				setPStar( conf.getMinEnergyLB() );
			}

			startTime = System.currentTimeMillis();

			// start conformation queue
			confs.start();

			if( waitUntilCapacity )
				confs.waitUntilCapacity();

		}  catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}

	}


	protected void iterate() throws Exception {

		for( ServerInterface serverInterface : serverInterfaces ) {

			if( serverInterface.isProcessed() ) {

				// getProcessedConfs sets server to the ready state
				accumulate(serverInterface.getProcessedConfs());

				if( eAppx == EApproxReached.FALSE ) {
					// we are finished, so other threads can terminate
					break;
				}

			}

			if( serverInterface.isReady() ) {
				// server is idle and awaiting new conformations

				// lock conformation queue
				synchronized( confs.qLock ) {

					int request = unProcessedConfs.size();
					int granted = 0;

					if( (granted = canSatisfy(request)) == 0 )
						break;

					// reduce the size of buf to match
					while( unProcessedConfs.size() > granted ) {
						unProcessedConfs.remove(unProcessedConfs.size()-1);
					}

					for( int i = 0; i < unProcessedConfs.size(); ++i ) 
						unProcessedConfs.set(i, confs.deQueue());

					minimizingConfs = minimizingConfs.add( BigInteger.valueOf(unProcessedConfs.size()) );

					if( confs.getState() == Thread.State.WAITING ) confs.qLock.notify();
				}

				// send in unprocessed confs
				serverInterface.setUnprocessedConfs( unProcessedConfs );
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

				cleanUpSlaves( serverInterfaces );
			}

		}  catch (Exception e) {
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

			cleanUpSlaves( serverInterfaces );

		}  catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	public void cleanUpSlaves( ArrayList<ServerInterface> serverInterfaces ) throws RemoteException {
		try {

			for( ServerInterface serverInterface : serverInterfaces ) {

				//while( eAppx != EApproxReached.TRUE && !serverInterface.isReady() ) {
				if( eAppx != EApproxReached.TRUE && !serverInterface.isReady() ) {
					// if server is not ready, then it is processing conformations
					while( !serverInterface.isProcessed() ) Thread.sleep(sleepInterval);

					accumulate(serverInterface.getProcessedConfs());
				}

				serverInterface.softTerminate();
			}
		} catch( Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}

}