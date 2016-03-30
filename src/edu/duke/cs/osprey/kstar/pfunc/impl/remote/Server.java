package edu.duke.cs.osprey.kstar.pfunc.impl.remote;

import java.net.MalformedURLException;
import java.net.UnknownHostException;
import java.rmi.AlreadyBoundException;
import java.rmi.RemoteException;
import java.rmi.registry.LocateRegistry;
import java.rmi.registry.Registry;

public class Server {

	public static void main(String args[]) throws RemoteException, 
	MalformedURLException, AlreadyBoundException, InterruptedException, UnknownHostException {

		ServerImplementation serverImplementation = new ServerImplementation();
		Registry registry = LocateRegistry.createRegistry(Constants.RMI_PORT);
		registry.bind(Constants.RMI_ID, serverImplementation);

		System.out.println("Server:\t" + Constants.RMI_ID + " is started on " 
				+  java.net.InetAddress.getLocalHost().getHostName());
	}
}
