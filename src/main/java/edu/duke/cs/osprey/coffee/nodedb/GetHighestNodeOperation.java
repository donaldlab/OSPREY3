package edu.duke.cs.osprey.coffee.nodedb;

import com.hazelcast.nio.ObjectDataInput;
import com.hazelcast.nio.ObjectDataOutput;
import com.hazelcast.spi.impl.operationservice.Operation;

import java.io.IOException;


/**
 * Returns the highest-scoring node from a cluster member.
 */
public class GetHighestNodeOperation extends Operation {

	private int statei;

	private NodeIndex.Node node = null;

	@SuppressWarnings("unused") // used by hazelcast
	public GetHighestNodeOperation() {
		this.statei = -1;
	}

	public GetHighestNodeOperation(int statei) {
		this.statei = statei;
	}

	@Override
	public final boolean returnsResponse() {
		return true;
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
		node = nodedb.removeHighestLocal(statei);
	}

	@Override
	public Object getResponse() {
		return node;
	}
}
