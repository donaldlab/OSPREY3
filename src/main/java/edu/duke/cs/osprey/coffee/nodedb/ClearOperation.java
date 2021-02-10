package edu.duke.cs.osprey.coffee.nodedb;

import com.hazelcast.nio.ObjectDataInput;
import com.hazelcast.nio.ObjectDataOutput;
import com.hazelcast.spi.impl.operationservice.Operation;

import java.io.IOException;


/**
 * Deletes all nodes for a state
 */
public class ClearOperation extends Operation {

	private int statei;

	@SuppressWarnings("unused") // used by hazelcast
	public ClearOperation() {
		this.statei = -1;
	}

	public ClearOperation(int statei) {
		this.statei = statei;
	}

	@Override
	public final boolean returnsResponse() {
		return false;
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
	public String getServiceName() {
		return NodeDB.ServiceName;
	}

	@Override
	public final void run() {
		NodeDB nodedb = getService();
		nodedb.clearLocal(statei);
	}
}
