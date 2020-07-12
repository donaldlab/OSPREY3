package edu.duke.cs.osprey.coffee.nodedb;

import com.hazelcast.nio.ObjectDataInput;
import com.hazelcast.nio.ObjectDataOutput;
import com.hazelcast.spi.impl.operationservice.Operation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


/**
 * Returns the highest-scoring node from a cluster member.
 */
public class GetHighestNodesOperation extends Operation {

	private int statei;
	private int count;

	private List<NodeIndex.Node> nodes = null;

	@SuppressWarnings("unused") // used by hazelcast
	public GetHighestNodesOperation() {
		this.statei = -1;
		this.count = -1;
	}

	public GetHighestNodesOperation(int statei, int count) {
		this.statei = statei;
		this.count = count;
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
		out.writeInt(count);
	}

	@Override
	protected void readInternal(ObjectDataInput in)
	throws IOException {
		super.readInternal(in);
		statei = in.readInt();
		count = in.readInt();
	}

	@Override
	public String getServiceName() {
		return NodeDB.ServiceName;
	}

	@Override
	public final void run() {
		NodeDB nodedb = getService();
		nodes = new ArrayList<>(count);
		nodedb.removeHighestLocal(statei, count, nodes);
	}

	@Override
	public Object getResponse() {
		return nodes;
	}
}
