package edu.duke.cs.osprey.coffee.directions;

import com.hazelcast.nio.ObjectDataInput;
import com.hazelcast.nio.ObjectDataOutput;
import com.hazelcast.spi.impl.operationservice.Operation;
import edu.duke.cs.osprey.confspace.Sequence;

import java.io.IOException;
import java.util.Collection;


public class FinishedSequencesOperation extends Operation {

	private int sequencedStatei;
	private int[][] sequences;

	@SuppressWarnings("unused") // used by hazelcast
	public FinishedSequencesOperation() {
		sequencedStatei = -1;
		sequences = null;
	}

	public FinishedSequencesOperation(int sequencedStatei, Collection<Sequence> sequences) {

		if (sequences.isEmpty()) {
			throw new IllegalArgumentException("need sequences to finish");
		}

		this.sequencedStatei = sequencedStatei;
		this.sequences = new int[sequences.size()][];
		int i = 0;
		for (var seq : sequences) {
			this.sequences[i++] = seq.rtIndices;
		}
	}

	@Override
	public final boolean returnsResponse() {
		return false;
	}

	@Override
	public String getServiceName() {
		return Directions.ServiceName;
	}

	@Override
	protected void writeInternal(ObjectDataOutput out)
	throws IOException {
		super.writeInternal(out);

		int numPos = sequences[0].length;

		out.writeInt(sequencedStatei);
		out.writeInt(sequences.length);
		out.writeInt(numPos);
		for (var seq : sequences) {
			assert (seq.length == numPos);
			for (var i : seq) {
				out.writeInt(i);
			}
		}
	}

	@Override
	protected void readInternal(ObjectDataInput in)
	throws IOException {
		super.readInternal(in);

		sequencedStatei = in.readInt();
		int numSeqs = in.readInt();
		int numPos = in.readInt();
		sequences = new int[numSeqs][numPos];
		for (int i=0; i<numSeqs; i++) {
			for (int posi=0; posi<numPos; posi++) {
				sequences[i][posi] = in.readInt();
			}
		}
	}

	@Override
	public final void run() {
		Directions directions = getService();
		directions.receiveFinishedSequences(sequencedStatei, sequences);
	}
}
