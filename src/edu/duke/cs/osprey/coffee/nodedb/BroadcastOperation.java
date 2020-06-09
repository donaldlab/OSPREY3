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

	private long[] freeSpaces;
	private BigExp[] maxScores;

	@SuppressWarnings("unused") // used by hazelcast
	public BroadcastOperation() {
		this.freeSpaces = null;
		this.maxScores = null;
	}

	public BroadcastOperation(long[] freeSpaces, BigExp[] maxScores) {
		this.freeSpaces = freeSpaces;
		this.maxScores = maxScores;
	}

	@Override
	public final boolean returnsResponse() {
		return false;
	}

	@Override
	protected void writeInternal(ObjectDataOutput out)
	throws IOException {
		super.writeInternal(out);

		assert (freeSpaces.length == maxScores.length);
		int numStates = freeSpaces.length;
		out.writeInt(numStates);

		for (long freeSpace : freeSpaces) {
			out.writeLong(freeSpace);
		}

		for (var maxScore : maxScores) {
			if (maxScore != null) {
				out.writeBoolean(true);
				out.writeDouble(maxScore.fp);
				out.writeInt(maxScore.exp);
			} else {
				out.writeBoolean(false);
			}
		}
	}

	@Override
	protected void readInternal(ObjectDataInput in)
	throws IOException {
		super.readInternal(in);

		int numStates = in.readInt();
		freeSpaces = new long[numStates];
		maxScores = new BigExp[numStates];

		for (int statei=0; statei<maxScores.length; statei++) {
			freeSpaces[statei] = in.readLong();
		}

		for (int statei=0; statei<maxScores.length; statei++) {
			maxScores[statei] = null;
			if (in.readBoolean()) {
				double fp = in.readDouble();
				int exp = in.readInt();
				maxScores[statei] = new BigExp(fp, exp);
			}
		}
	}

	@Override
	public String getServiceName() {
		return NodeDB.ServiceName;
	}

	@Override
	public final void run() {
		NodeDB nodedb = getService();
		nodedb.receiveBroadcast(getCallerAddress(), freeSpaces, maxScores);
	}
}
