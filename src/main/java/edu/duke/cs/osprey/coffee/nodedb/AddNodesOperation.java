package edu.duke.cs.osprey.coffee.nodedb;

import com.hazelcast.nio.ObjectDataInput;
import com.hazelcast.nio.ObjectDataOutput;
import com.hazelcast.spi.impl.operationservice.Operation;
import edu.duke.cs.osprey.coffee.Serializers;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


/**
 * Adds nodes to the cluster member.
 */
public class AddNodesOperation extends Operation {

	private int statei;
	private List<NodeIndex.Node> nodes;

	@SuppressWarnings("unused") // used by hazelcast
	public AddNodesOperation() {
		this.statei = -1;
		this.nodes = null;
	}

	public AddNodesOperation(int statei, List<NodeIndex.Node> nodes) {
		this.statei = statei;
		this.nodes = nodes;
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
		out.writeInt(statei);
		out.writeInt(nodes.size());
		for (var node : nodes) {
			Serializers.hazelcastNodeSerializeDeNovo(out, nodedb.confSpace, node);
		}
	}

	@Override
	protected void readInternal(ObjectDataInput in)
	throws IOException {
		super.readInternal(in);

		// NOTE: can't getService() here  ;_;

		statei = in.readInt();
		int size = in.readInt();
		nodes = new ArrayList<>(size);
		for (int i=0; i<size; i++) {
			nodes.add(Serializers.hazelcastNodeDeserializeDeNovo(in));
		}
	}

	@Override
	public String getServiceName() {
		return NodeDB.ServiceName;
	}

	@Override
	public final void run() {
		NodeDB nodedb = getService();
		nodedb.addLocal(statei, nodes);
	}
}
