package edu.duke.cs.osprey.coffee.nodedb;

import com.hazelcast.nio.ObjectDataInput;
import com.hazelcast.nio.ObjectDataOutput;
import com.hazelcast.spi.impl.operationservice.Operation;
import edu.duke.cs.osprey.coffee.Serializers;

import java.io.IOException;


/**
 * Adds a node to the cluster member.
 */
public class AddNodeOperation extends Operation {

	private NodeIndex.Node node;

	@SuppressWarnings("unused") // used by hazelcast
	public AddNodeOperation() {
		this.node = null;
	}

	public AddNodeOperation(NodeIndex.Node node) {
		this.node = node;
	}

	@Override
	public final boolean returnsResponse() {
		return false;
	}

	@Override
	protected void writeInternal(ObjectDataOutput out)
	throws IOException {
		super.writeInternal(out);

		NodeDB nodedb = getService();
		Serializers.hazelcastNodeSerializeDeNovo(out, nodedb.confSpace, node);
	}

	@Override
	protected void readInternal(ObjectDataInput in)
	throws IOException {
		super.readInternal(in);

		// NOTE: can't getService() here  ;_;

		node = Serializers.hazelcastNodeDeserializeDeNovo(in);
	}

	@Override
	public String getServiceName() {
		return NodeDB.ServiceName;
	}

	@Override
	public final void run() {
		NodeDB nodedb = getService();
		nodedb.addLocal(node);
	}
}
