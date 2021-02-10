package edu.duke.cs.osprey.coffee.zmat;

import com.hazelcast.nio.ObjectDataInput;
import com.hazelcast.nio.ObjectDataOutput;
import com.hazelcast.spi.impl.operationservice.Operation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


public class TuplesOperation extends Operation {

	public List<Tuple> tuples;

	@SuppressWarnings("unused") // used by hazelcast
	public TuplesOperation() {
		tuples = null;
	}

	public TuplesOperation(List<Tuple> tuples) {
		this.tuples = tuples;
	}

	@Override
	public final boolean returnsResponse() {
		return false;
	}

	@Override
	public String getServiceName() {
		return ClusterZMatrix.ServiceName;
	}

	@Override
	protected void writeInternal(ObjectDataOutput out)
	throws IOException {
		super.writeInternal(out);

		// all the tuples should be the same type
		if (tuples.stream().anyMatch(tuple -> tuple.type() != tuples.get(0).type())) {
			throw new IllegalArgumentException("heterogeneous tuple batches not supported");
		}

		out.writeInt(tuples.size());
		out.writeInt(tuples.get(0).type());

		for (var tuple : tuples) {
			tuple.write(out);
		}
	}

	@Override
	protected void readInternal(ObjectDataInput in)
	throws IOException {
		super.readInternal(in);

		int size = in.readInt();
		int type = in.readInt();

		tuples = new ArrayList<>(size);
		for (int i=0; i<size; i++) {
			switch (type) {
				case Single.Type -> tuples.add(Single.read(in));
				case Pair.Type -> tuples.add(Pair.read(in));
				case Triple.Type -> tuples.add(Triple.read(in));
				default -> throw new IllegalArgumentException("unrecognized tuple type: " + type);
			}
		}
	}

	@Override
	public final void run() {
		ClusterZMatrix zmat = getService();
		zmat.writeTuples(tuples);
	}
}
