package edu.duke.cs.osprey.kstar.pfunc.impl.remote;

import java.rmi.Remote;
import java.rmi.RemoteException;
import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;

public interface ServerInterface extends Remote {

	public boolean envVarsSet() throws RemoteException;
	
	public void setEnvVars( int fibers, int threadsPerFiber, 
			int confsPerThread,
			String eMinMethod,
			EnergyFunctionGenerator curEFcnGenerator,
			GenericResidueTemplateLibrary resTemplates,
			boolean assignTemplatesToStruct,
			boolean deleteNonTemplateResidues,
			boolean useMPI,
			double DUNBRACK_PROBABILTY_CUTOFF,
			String dataDir ) throws RemoteException;
	
	public void initFibers( SearchProblem sp ) throws RemoteException;

	public void setUnprocessedConfs( ArrayList<KSConf> confs ) throws RemoteException;

	public ArrayList<KSConf> getProcessedConfs() throws RemoteException;

	public boolean isProcessed() throws RemoteException;

	public boolean isReady() throws RemoteException;
	
	public void setEpsilon( boolean epsilon ) throws RemoteException;

	public void softTerminate() throws RemoteException;

	public void hardTerminate() throws RemoteException;

}
