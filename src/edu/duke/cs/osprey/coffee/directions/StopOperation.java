package edu.duke.cs.osprey.coffee.directions;

import com.hazelcast.spi.impl.operationservice.Operation;


public class StopOperation extends Operation {

	@Override
	public final boolean returnsResponse() {
		return false;
	}

	@Override
	public String getServiceName() {
		return Directions.ServiceName;
	}

	@Override
	public final void run() {
		Directions directions = getService();
		directions.receiveStop();
	}
}
