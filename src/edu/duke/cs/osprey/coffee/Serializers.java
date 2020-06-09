package edu.duke.cs.osprey.coffee;

import com.hazelcast.nio.ObjectDataInput;
import com.hazelcast.nio.ObjectDataOutput;
import com.hazelcast.nio.serialization.StreamSerializer;
import edu.duke.cs.osprey.coffee.nodedb.FixedIndex;
import edu.duke.cs.osprey.coffee.nodedb.NodeIndex;
import edu.duke.cs.osprey.coffee.seqdb.SaveOperation;
import edu.duke.cs.osprey.coffee.seqdb.SeqInfo;
import edu.duke.cs.osprey.confspace.ConfSpaceIteration;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.IntEncoding;
import edu.duke.cs.osprey.tools.MapDBTools;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import org.jetbrains.annotations.NotNull;
import org.mapdb.DataInput2;
import org.mapdb.DataOutput2;

import java.io.IOException;
import java.math.BigDecimal;
import java.nio.ByteBuffer;
import java.util.Comparator;
import java.util.stream.IntStream;


public class Serializers {

	public static final MapDBTools.BigExpSerializer BigExp = new MapDBTools.BigExpSerializer();

	private static int getMaxConf(ConfSpaceIteration confSpace) {
		return IntStream.range(0, confSpace.numPos())
			.map(posi -> confSpace.numConf(posi) - 1)
			.max()
			.orElseThrow();
	}

	private static IntEncoding getConfEncoding(ConfSpaceIteration confSpace) {
		return IntEncoding.get(writeConf(getMaxConf(confSpace)));
	}

	private static int writeConf(int conf) {
		// add 1 to the max value to map the range [-1,N) to [0,N]
		// since -1 is used for unassigned values
		return conf + 1;
	}

	private static int readConf(int conf) {
		// undo the effects of writeConf
		return conf - 1;
	}

	private static int getMaxSeq(SeqSpace seqSpace) {
		return seqSpace.positions.stream()
			.mapToInt(pos -> pos.resTypes.size() - 1)
			.max()
			.orElseThrow();
	}

	private static IntEncoding getSeqEncoding(SeqSpace seqSpace) {
		return IntEncoding.get(writeSeq(getMaxSeq(seqSpace)));
	}

	private static int writeSeq(int seq) {
		// add 1 to the max value to map the range [-1,N) to [0,N]
		// since -1 is used for unassigned values
		return seq + 1;
	}

	private static int readSeq(int seq) {
		// undo the effects of writeSeq
		return seq - 1;
	}

	public static MapDBTools.IntArraySerializer mapdbConf(ConfSpaceIteration confSpace) {
		return new MapDBTools.IntArraySerializer(getMaxConf(confSpace), confSpace.numPos());
	}

	public static MapDBTools.IntArraySerializer mapdbSeq(SeqSpace seqSpace) {
		return new MapDBTools.IntArraySerializer(
			seqSpace.positions.stream()
				.flatMap(pos -> pos.resTypes.stream())
				.mapToInt(resType -> resType.index)
				.max()
				.orElseThrow(),
			seqSpace.positions.size()
		);
	}

	public static FixedIndex.Serializer<NodeIndex.Node> indexNode(MultiStateConfSpace.State state) {

		var confEncoding = getConfEncoding(state.confSpace);

		return new FixedIndex.Serializer<>() {

			@Override
			public int bytes() {
				return 8+4 // BigExp
					+ state.confSpace.numPos()*confEncoding.numBytes; // conf
			}

			@Override
			public void serialize(ByteBuffer out, NodeIndex.Node node) {
				out.putDouble(node.score.fp);
				out.putInt(node.score.exp);
				for (int i=0; i<state.confSpace.numPos(); i++) {
					confEncoding.write(out, writeConf(node.conf[i]));
				}
			}

			@Override
			public NodeIndex.Node deserialize(ByteBuffer in) {
				double fp = in.getDouble();
				int exp = in.getInt();
				BigExp score = new BigExp(fp, exp);
				int[] conf = new int[state.confSpace.numPos()];
				for (int i=0; i<conf.length; i++) {
					conf[i] = readConf(confEncoding.read(in));
				}
				return new NodeIndex.Node(state.index, conf, score);
			}
		};
	}

	protected static abstract class AbstractStreamSerializer<T> implements StreamSerializer<T> {

		public final int typeId;

		protected AbstractStreamSerializer(int typeId) {
			this.typeId = typeId;
		}

		@Override
		public int getTypeId() {
			return typeId;
		}

		@Override
		public void destroy() {
			// nothing needed here
		}
	}

	// don't collide with hazelcasts's serialization IDs in SerializationConstants
	// any positive ID seems to be safe
	public static final int IdNode = 1001;

	public static StreamSerializer<NodeIndex.Node> hazelcastNode(MultiStateConfSpace confSpace) {

		var stateEncoding = IntEncoding.get(confSpace.states.size() - 1);

		// find the biggest conf encoding needed for all of the states
		var confEncoding = confSpace.states.stream()
			.map(state -> getConfEncoding(state.confSpace))
			.max(Comparator.comparing(encoding -> encoding.numBytes))
			.orElseThrow();

		return new AbstractStreamSerializer<>(IdNode) {

			@Override
			public void write(ObjectDataOutput out, NodeIndex.Node node)
			throws IOException {
				stateEncoding.write(out, node.statei);
				out.writeDouble(node.score.fp);
				out.writeInt(node.score.exp);
				int numPos = confSpace.states.get(node.statei).confSpace.numPos();
				for (int posi=0; posi<numPos; posi++) {
					confEncoding.write(out, writeConf(node.conf[posi]));
				}
			}

			@Override
			public NodeIndex.Node read(ObjectDataInput in)
			throws IOException {
				int statei = stateEncoding.read(in);
				double fp = in.readDouble();
				int exp = in.readInt();
				BigExp score = new BigExp(fp, exp);
				int numPos = confSpace.states.get(statei).confSpace.numPos();
				int[] conf = new int[numPos];
				for (int posi=0; posi<numPos; posi++) {
					conf[posi] = readConf(confEncoding.read(in));
				}
				return new NodeIndex.Node(statei, conf, score);
			}
		};
	}

	public static void hazelcastNodeSerializeDeNovo(ObjectDataOutput out, MultiStateConfSpace confSpace, NodeIndex.Node node)
	throws IOException {

		// need to serialize so the other end can deserialize without any context
		// so we can't use any fancy compression techniques based on knowing details of the conf space =(

		out.writeInt(node.statei);
		out.writeDouble(node.score.fp);
		out.writeInt(node.score.exp);
		int numPos = confSpace.states.get(node.statei).confSpace.numPos();
		out.writeInt(numPos);
		for (int posi=0; posi<numPos; posi++) {
			out.writeInt(writeConf(node.conf[posi]));
		}
	}

	public static NodeIndex.Node hazelcastNodeDeserializeDeNovo(ObjectDataInput in)
	throws IOException {

		int statei = in.readInt();
		double fp = in.readDouble();
		int exp = in.readInt();
		BigExp score = new BigExp(fp, exp);
		int numPos = in.readInt();
		int[] conf = new int[numPos];
		for (int posi=0; posi<numPos; posi++) {
			conf[posi] = readConf(in.readInt());
		}

		return new NodeIndex.Node(statei, conf, score);
	}


	public static MapDBTools.SimpleSerializer<SeqInfo> mapdbSeqInfo(MultiStateConfSpace confSpace) {
		return new MapDBTools.SimpleSerializer<>() {

			private final MapDBTools.BigDecimalBoundsSerializer s = new MapDBTools.BigDecimalBoundsSerializer();

			@Override
			public void serialize(@NotNull DataOutput2 out, @NotNull SeqInfo data)
				throws IOException {
				for (int statei=0; statei<confSpace.sequencedStates.size(); statei++) {
					s.serialize(out, data.zSumBounds[statei]);
				}
			}

			@Override
			public SeqInfo deserialize(@NotNull DataInput2 in, int available)
				throws IOException {
				SeqInfo data = new SeqInfo(confSpace);
				for (int statei=0; statei<confSpace.sequencedStates.size(); statei++) {
					data.zSumBounds[statei] = s.deserialize(in, available);
				}
				return data;
			}

			@Override
			public int compare(SeqInfo a, SeqInfo b) {
				throw new UnsupportedOperationException();
			}
		};
	}

	public static void hazelcastSeqDBSerializeDeNovo(ObjectDataOutput out, MultiStateConfSpace confSpace, SaveOperation op)
	throws IOException {

		// need to serialize so the other end can deserialize without any context

		// write sequence space metadata
		int numSequencedStates = confSpace.sequencedStates.size();
		out.writeInt(numSequencedStates);
		int numPos = confSpace.seqSpace.positions.size();
		out.writeInt(numPos);
		var seqEncoding = getSeqEncoding(confSpace.seqSpace);
		out.writeByte(seqEncoding.ordinal());
		var stateEncoding = IntEncoding.get(confSpace.sequencedStates.size() - 1);
		out.writeByte(stateEncoding.ordinal());

		// write the sequenced sums
		out.writeInt(op.sequencedSums.length);
		for (int sumi=0; sumi<op.sequencedSums.length; sumi++) {
			var sum = op.sequencedSums[sumi];
			for (int posi=0; posi<numPos; posi++) {
				seqEncoding.write(out, writeSeq(sum.seq[posi]));
			}
			for (int statei=0; statei<numSequencedStates; statei++) {
				out.writeObject(sum.boundsByState[statei].lower);
				out.writeObject(sum.boundsByState[statei].upper);
			}
		}

		// write the unsequenced sums
		out.writeInt(op.unsequencedSums.length);
		for (int sumi=0; sumi<op.unsequencedSums.length; sumi++) {
			var sum = op.unsequencedSums[sumi];
			stateEncoding.write(out, sum.unsequencedIndex);
			out.writeObject(sum.bounds.lower);
			out.writeObject(sum.bounds.upper);
		}
	}

	public static void hazelcastSeqDBDeserializeDeNovo(ObjectDataInput in, SaveOperation op)
	throws IOException {

		// read the sequence space metadata
		int numSequencedStates = in.readInt();
		int numPos = in.readInt();
		IntEncoding seqEncoding = IntEncoding.values()[in.readByte()];
		IntEncoding stateEncoding = IntEncoding.values()[in.readByte()];

		// read the sequenced sums
		op.sequencedSums = new SaveOperation.SequencedSum[in.readInt()];
		for (int sumi=0; sumi<op.sequencedSums.length; sumi++) {
			int[] seq = new int[numPos];
			for (int posi=0; posi<numPos; posi++) {
				seq[posi] = readSeq(seqEncoding.read(in));
			}
			var boundsByState = new BigDecimalBounds[numSequencedStates];
			for (int statei=0; statei<numSequencedStates; statei++) {
				BigDecimal lower = in.readObject();
				BigDecimal upper = in.readObject();
				boundsByState[statei] = new BigDecimalBounds(lower, upper);
			}
			op.sequencedSums[sumi] = new SaveOperation.SequencedSum(seq, boundsByState);
		}

		// read the unsequenced sums
		op.unsequencedSums = new SaveOperation.UnsequencedSum[in.readInt()];
		for (int sumi=0; sumi<op.unsequencedSums.length; sumi++) {
			int sequencedIndex = stateEncoding.read(in);
			BigDecimal lower = in.readObject();
			BigDecimal upper = in.readObject();
			op.unsequencedSums[sumi] = new SaveOperation.UnsequencedSum(sequencedIndex, new BigDecimalBounds(lower, upper));
		}
	}
}
