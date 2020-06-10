package edu.duke.cs.osprey.coffee.commands;

import com.hazelcast.nio.ObjectDataInput;
import com.hazelcast.nio.ObjectDataOutput;
import com.hazelcast.spi.impl.operationservice.Operation;

import java.io.IOException;


public class FocusOperation extends Operation {

	private int statei;

	public FocusOperation() {
		statei = -1;
	}

	public FocusOperation(int statei) {
		this.statei = statei;
	}

	@Override
	public final boolean returnsResponse() {
		return false;
	}

	@Override
	public String getServiceName() {
		return Commands.ServiceName;
	}

	@Override
	protected void writeInternal(ObjectDataOutput out)
	throws IOException {
		super.writeInternal(out);

		out.writeInt(statei);
	}

	@Override
	protected void readInternal(ObjectDataInput in)
	throws IOException {
		super.readInternal(in);

		statei = in.readInt();
	}

	@Override
	public final void run() {
		Commands commands = getService();
		commands.receiveFocus(statei);
	}
}
