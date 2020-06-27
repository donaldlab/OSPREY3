package edu.duke.cs.osprey.coffee.seqdb;

import com.hazelcast.nio.ObjectDataInput;
import com.hazelcast.nio.ObjectDataOutput;
import com.hazelcast.spi.impl.operationservice.Operation;
import edu.duke.cs.osprey.coffee.Serializers;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;

import java.io.IOException;


public class SaveOperation extends Operation {

	public static class SequencedSum {

		public final int[] seq;
		public final StateZ[] statezs;

		public SequencedSum(int[] seq, StateZ[] statezs) {
			this.seq = seq;
			this.statezs = statezs;
		}
	}

	public static class UnsequencedSum {

		public final StateZ statez;

		public UnsequencedSum(StateZ statez) {
			this.statez = statez;
		}
	}

	public int numSequencedStates;
	public SequencedSum[] sequencedSums;
	public UnsequencedSum[] unsequencedSums;

	public SaveOperation() {
		numSequencedStates = 0;
		sequencedSums = null;
		unsequencedSums = null;
	}

	public SaveOperation(MultiStateConfSpace confSpace, Batch batch) {

		numSequencedStates = confSpace.sequencedStates.size();

		sequencedSums = batch.sequencedSums.entrySet().stream()
			.map(entry -> {
				Sequence seq = entry.getKey();
				SeqInfo info = entry.getValue();
				assert (info.statezs.length == numSequencedStates);
				return new SequencedSum(seq.rtIndices, info.statezs);
			})
			.toArray(SequencedSum[]::new);

		unsequencedSums = batch.unsequencedSums.values().stream()
			.map(statez -> new UnsequencedSum(statez))
			.toArray(UnsequencedSum[]::new);
	}

	@Override
	public final boolean returnsResponse() {
		return false;
	}

	@Override
	public String getServiceName() {
		return SeqDB.ServiceName;
	}

	@Override
	protected void writeInternal(ObjectDataOutput out)
	throws IOException {
		super.writeInternal(out);

		SeqDB seqdb = getService();
		Serializers.hazelcastSeqDBSerializeDeNovo(out, seqdb.confSpace, this);
	}

	@Override
	protected void readInternal(ObjectDataInput in)
	throws IOException {
		super.readInternal(in);

		Serializers.hazelcastSeqDBDeserializeDeNovo(in, this);
	}

	@Override
	public final void run() {
		SeqDB seqdb = getService();
		seqdb.commitBatch(this);
	}
}
