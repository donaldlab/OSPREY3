package edu.duke.cs.osprey.coffee.nodedb;

import com.hazelcast.nio.ObjectDataInput;
import com.hazelcast.nio.ObjectDataOutput;
import com.hazelcast.spi.impl.operationservice.Operation;
import edu.duke.cs.osprey.tools.BigExp;

import java.io.IOException;


/**
 * Informs other cluster members about a member's current state.
 */
public class BroadcastOperation extends Operation {

	private NodeIndices.BroadcastInfo nodesInfo;

	// TODO: serialize node performance

	@SuppressWarnings("unused") // used by hazelcast
	public BroadcastOperation() {
		nodesInfo = null;
	}

	public BroadcastOperation(NodeIndices.BroadcastInfo nodesInfo, NodePerformance nodePerformance) {
		this.nodesInfo = nodesInfo;
	}

	@Override
	public final boolean returnsResponse() {
		return false;
	}

	@Override
	protected void writeInternal(ObjectDataOutput out)
	throws IOException {
		super.writeInternal(out);

		int n = nodesInfo.size();
		out.writeInt(n);

		for (var freeSpace : nodesInfo.freeSpaces) {
			out.writeLong(freeSpace);
		}

		for (var maxScore : nodesInfo.maxScores) {
			if (maxScore != null) {
				out.writeBoolean(true);
				out.writeDouble(maxScore.fp);
				out.writeInt(maxScore.exp);
			} else {
				out.writeBoolean(false);
			}
		}

		out.writeLong(nodesInfo.usedBytes);
		out.writeLong(nodesInfo.totalBytes);
	}

	@Override
	protected void readInternal(ObjectDataInput in)
	throws IOException {
		super.readInternal(in);

		int n = in.readInt();
		nodesInfo = new NodeIndices.BroadcastInfo(n);

		for (int i=0; i<n; i++) {
			nodesInfo.freeSpaces[i] = in.readLong();
		}

		for (int i=0; i<n; i++) {
			nodesInfo.maxScores[i] = null;
			if (in.readBoolean()) {
				double fp = in.readDouble();
				int exp = in.readInt();
				nodesInfo.maxScores[i] = new BigExp(fp, exp);
			}
		}

		nodesInfo.usedBytes = in.readLong();
		nodesInfo.totalBytes = in.readLong();
	}

	@Override
	public String getServiceName() {
		return NodeDB.ServiceName;
	}

	@Override
	public final void run() {
		NodeDB nodedb = getService();
		nodedb.receiveBroadcast(getCallerAddress(), nodesInfo);
	}
}
