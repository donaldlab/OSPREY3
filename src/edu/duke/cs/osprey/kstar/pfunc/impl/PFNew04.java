package edu.duke.cs.osprey.kstar.pfunc.impl;

import java.io.Serializable;
import java.math.BigInteger;
import java.rmi.RemoteException;
import java.rmi.registry.LocateRegistry;
import java.rmi.registry.Registry;
import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.kstar.KSConfQ;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.impl.remote.Constants;
import edu.duke.cs.osprey.kstar.pfunc.impl.remote.ServerInterface;
import edu.duke.cs.osprey.kstar.KSConf;

@SuppressWarnings("serial")
public class PFNew04 extends PFNew03 implements Serializable {

	ArrayList<ServerInterface> serverInterfaces = new ArrayList<>();
	ArrayList<KSConf> unProcessedConfs = new ArrayList<>();

	public PFNew04() {
		super();
	}

	public PFNew04( int strand, ArrayList<String> sequence, 
			ArrayList<Integer> absolutePos, 
			String checkPointPath, String reducedSPName, 
			ConfigFileParser cfp, SearchProblem panSP ) {

		super( strand, sequence, absolutePos, checkPointPath, reducedSPName, cfp, panSP );
	}

	public void start() {

		try {

			setRunState(RunState.STARTED);
			
			// set pstar
			initTradPStar();

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
				serverInterface.initFibers(reducedSP);
			}

			unProcessedConfs.clear();
			for( int i = 0; i < fibers * confsPerThread; ++i ) unProcessedConfs.add(null);
			unProcessedConfs.trimToSize();

			confsQ = new KSConfQ( this, unProcessedConfs.size(), partialQLB );

			// start conformation queue
			confsQ.start();

		}  catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}

		startTime = System.currentTimeMillis();
	}


	protected void iterate() throws Exception {

		for( ServerInterface serverInterface : serverInterfaces ) {

			if( serverInterface.isProcessed() ) {

				// getProcessedConfs sets server to the ready state
				accumulate(serverInterface.getProcessedConfs(), true);

				if( eAppx == EApproxReached.FALSE ) {
					// we are finished, so other threads can terminate
					break;
				}

			}

			if( serverInterface.isReady() ) {
				// server is idle and awaiting new conformations

				// lock conformation queue
				synchronized( confsQ.lock ) {

					int request = unProcessedConfs.size();
					int granted = 0;

					if( (granted = canSatisfy(request)) == 0 )
						break;

					// reduce the size of buf to match
					while( unProcessedConfs.size() > granted ) {
						unProcessedConfs.remove(unProcessedConfs.size()-1);
					}

					for( int i = 0; i < unProcessedConfs.size(); ++i ) 
						unProcessedConfs.set(i, confsQ.deQueue());

					minimizingConfs = minimizingConfs.add( BigInteger.valueOf(unProcessedConfs.size()) );

					if( confsQ.getState() == Thread.State.WAITING ) confsQ.lock.notify();
				}

				// send in unprocessed confs
				serverInterface.setUnprocessedConfs( unProcessedConfs );
			}
		}
	}


	protected void computeSlice() {

		try {

			iterate();

			if( eAppx == EApproxReached.FALSE ) Thread.sleep(sleepInterval);

			else {
				confsQ.cleanUp(true);

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

			confsQ.cleanUp(true);

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

					accumulate(serverInterface.getProcessedConfs(), true);
				}

				serverInterface.softTerminate();
			}
		} catch( Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	public String getImpl() {
		return "new04";
	}

}