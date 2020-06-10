package edu.duke.cs.osprey.coffee.commands;

import com.hazelcast.spi.impl.operationservice.Operation;


public class StopOperation extends Operation {

	@Override
	public final boolean returnsResponse() {
		return false;
	}

	@Override
	public String getServiceName() {
		return Commands.ServiceName;
	}

	@Override
	public final void run() {
		Commands commands = getService();
		commands.receiveStop();
	}
}
