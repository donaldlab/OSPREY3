package edu.duke.cs.osprey.kstar.remote;

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
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.PF1NMTPCPMCache;
import edu.duke.cs.osprey.kstar.PFAbstract;
import edu.duke.cs.osprey.pruning.PruningControl;

public class PFMNPCPMCache extends PF1NMTPCPMCache {

	ArrayList<ServerInterface> serverInterfaces = new ArrayList<>();
	ArrayList<KSConf> unProcessedConfs = new ArrayList<>();


	public PFMNPCPMCache(ArrayList<String> sequence, ConfigFileParser cfp, 
			SearchProblem sp, PruningControl pc, DEEPerSettings dset, 
			ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones, 
			double EW_I0) {

		super( sequence, cfp, sp, pc, dset, moveableStrands, freeBBZones, EW_I0);

		// initialize environment variables for each server. this need only be done once
		try {

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

		}  catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}

	protected void start() {

		try {

			setRunState(RunState.STARTED);

			/*
			fibers = PFAbstract.getNumFibers();
			threadsPerFiber = PFAbstract.getNumThreads();
			confsPerThread = PFAbstract.getConfsThreadBuffer() * threadsPerFiber;

			// establish connection to minimizaiton servers
			serverInterfaces.clear();
			for( String server : serverList ) {

				System.out.println("Initiating connection to server: " + server);

				Registry registry = LocateRegistry.getRegistry( server, Constants.RMI_PORT );

				ServerInterface serverInterface = (ServerInterface)registry.lookup( Constants.RMI_ID );
				serverInterfaces.add( serverInterface );

				serverInterface.setEnvVars(fibers, threadsPerFiber, confsPerThread,
						PFAbstract.eMinMethod, 
						EnvironmentVars.curEFcnGenerator,
						EnvironmentVars.resTemplates,
						EnvironmentVars.assignTemplatesToStruct,
						EnvironmentVars.deleteNonTemplateResidues,
						EnvironmentVars.useMPI,
						EnvironmentVars.DUNBRACK_PROBABILTY_CUTOFF,
						EnvironmentVars.getDataDir());

				serverInterface.initFibers(sp);
			}
			 */

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

				setPStar( conf.getMinEnergyLowerBound() );
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
				if( ( eAppx = accumulate(serverInterface.getProcessedConfs()) ) != EApproxReached.FALSE ) {

					// we are finished, so other threads can terminate
					break;
				}

			}

			if( serverInterface.isReady() ) {
				// server is idle and awaiting new conformations

				// lock conformation queue
				synchronized( confs.qLock ) {
					// only process when there are enough confs ready to be processed
					if( confs.size() < unProcessedConfs.size() ) {

						if( !epsilonPossible(unProcessedConfs.size()) ) {

							break;
						}

						if( confs.size() < unProcessedConfs.size() ) confs.qLock.wait();
					}

					for( int i = 0; i < unProcessedConfs.size(); ++i ) unProcessedConfs.set(i, confs.deQueue());
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


	protected EApproxReached accumulate( ArrayList<KSConf> partialQConfs ) {

		for( KSConf conf : partialQConfs ) updateQStar( conf );

		// we need a current snapshot of qDagger, so we lock here
		synchronized( confs.qLock ) {

			minimizingConfs = minimizingConfs.subtract( BigInteger.valueOf(partialQConfs.size()) );

			confs.setQDagger( confs.getQDagger().subtract( computePartialQDagger(partialQConfs) ) );

			Et = confs.size() > 0 ? confs.get(confs.size()-1).getMinEnergyLowerBound() : partialQConfs.get(partialQConfs.size()-1).getMinEnergyLowerBound();

			// negative values of effective esilon are disallowed
			if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0) return EApproxReached.NOT_POSSIBLE;

			long currentTime = System.currentTimeMillis();

			System.out.println(partialQConfs.get(partialQConfs.size()-1).getMinEnergy() + "\t" + effectiveEpsilon + "\t" + 
					getNumMinimizedConfs() + "\t" + getNumUnMinimizedConfs() + "\t "+ ((currentTime-startTime)/1000));
		}

		return effectiveEpsilon > targetEpsilon ? EApproxReached.FALSE : EApproxReached.TRUE;
	}


	public void cleanUpSlaves( ArrayList<ServerInterface> serverInterfaces ) throws RemoteException {
		try {

			for( ServerInterface serverInterface : serverInterfaces ) {

				//while( eAppx != EApproxReached.TRUE && !serverInterface.isReady() ) {
				if( eAppx != EApproxReached.TRUE && !serverInterface.isReady() ) {
					// if server is not ready, then it is processing conformations
					while( !serverInterface.isProcessed() ) Thread.sleep(sleepInterval);

					eAppx = accumulate(serverInterface.getProcessedConfs());
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