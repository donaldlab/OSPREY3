package edu.duke.cs.osprey.coffee;

import com.hazelcast.nio.ObjectDataInput;
import com.hazelcast.nio.ObjectDataOutput;
import com.hazelcast.nio.serialization.StreamSerializer;
import edu.duke.cs.osprey.coffee.nodedb.FixedIndex;
import edu.duke.cs.osprey.coffee.nodedb.NodeIndex;
import edu.duke.cs.osprey.coffee.seqdb.SaveOperation;
import edu.duke.cs.osprey.coffee.seqdb.SeqInfo;
import edu.duke.cs.osprey.coffee.seqdb.StateZ;
import edu.duke.cs.osprey.confspace.ConfSearch;
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
import java.util.function.Consumer;
import java.util.function.Supplier;
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
				return state.confSpace.numPos()*confEncoding.numBytes // conf
					+ 8+4 // zSumUpper
					+ 8+4; // score
			}

			@Override
			public void serialize(ByteBuffer out, NodeIndex.Node node) {
				for (int i=0; i<state.confSpace.numPos(); i++) {
					confEncoding.write(out, writeConf(node.conf[i]));
				}
				out.putDouble(node.zSumUpper.fp);
				out.putInt(node.zSumUpper.exp);
				out.putDouble(node.score.fp);
				out.putInt(node.score.exp);
			}

			@Override
			public NodeIndex.Node deserialize(ByteBuffer in) {
				int[] conf = new int[state.confSpace.numPos()];
				for (int i=0; i<conf.length; i++) {
					conf[i] = readConf(confEncoding.read(in));
				}
				double fp = in.getDouble();
				int exp = in.getInt();
				BigExp zSumUpper = new BigExp(fp, exp);
				fp = in.getDouble();
				exp = in.getInt();
				BigExp score = new BigExp(fp, exp);
				return new NodeIndex.Node(state.index, conf, zSumUpper, score);
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
				int numPos = confSpace.states.get(node.statei).confSpace.numPos();
				for (int posi=0; posi<numPos; posi++) {
					confEncoding.write(out, writeConf(node.conf[posi]));
				}
				out.writeDouble(node.zSumUpper.fp);
				out.writeInt(node.zSumUpper.exp);
				out.writeDouble(node.score.fp);
				out.writeInt(node.score.exp);
			}

			@Override
			public NodeIndex.Node read(ObjectDataInput in)
			throws IOException {
				int statei = stateEncoding.read(in);
				int numPos = confSpace.states.get(statei).confSpace.numPos();
				int[] conf = new int[numPos];
				for (int posi=0; posi<numPos; posi++) {
					conf[posi] = readConf(confEncoding.read(in));
				}
				double fp = in.readDouble();
				int exp = in.readInt();
				BigExp zSumUpper = new BigExp(fp, exp);
				fp = in.readDouble();
				exp = in.readInt();
				BigExp score = new BigExp(fp, exp);
				return new NodeIndex.Node(statei, conf, zSumUpper, score);
			}
		};
	}

	public static void hazelcastNodeSerializeDeNovo(ObjectDataOutput out, MultiStateConfSpace confSpace, NodeIndex.Node node)
	throws IOException {

		// need to serialize so the other end can deserialize without any context
		// so we can't use any fancy compression techniques based on knowing details of the conf space =(

		out.writeInt(node.statei);
		int numPos = confSpace.states.get(node.statei).confSpace.numPos();
		out.writeInt(numPos);
		for (int posi=0; posi<numPos; posi++) {
			out.writeInt(writeConf(node.conf[posi]));
		}
		out.writeDouble(node.zSumUpper.fp);
		out.writeInt(node.zSumUpper.exp);
		out.writeDouble(node.score.fp);
		out.writeInt(node.score.exp);
	}

	public static NodeIndex.Node hazelcastNodeDeserializeDeNovo(ObjectDataInput in)
	throws IOException {

		int statei = in.readInt();
		int numPos = in.readInt();
		int[] conf = new int[numPos];
		for (int posi=0; posi<numPos; posi++) {
			conf[posi] = readConf(in.readInt());
		}
		double fp = in.readDouble();
		int exp = in.readInt();
		BigExp zSumUpper = new BigExp(fp, exp);
		fp = in.readDouble();
		exp = in.readInt();
		BigExp score = new BigExp(fp, exp);

		return new NodeIndex.Node(statei, conf, zSumUpper, score);
	}

	public static MapDBTools.SimpleSerializer<StateZ> mapdbStateZ(MultiStateConfSpace confSpace) {
		return new MapDBTools.SimpleSerializer<>() {

			private final MapDBTools.BigDecimalSerializer s = new MapDBTools.BigDecimalSerializer();
			private final IntEncoding[] confEncodings = confSpace.states.stream()
				.map(state -> getConfEncoding(state.confSpace))
				.toArray(IntEncoding[]::new);

			@Override
			public void serialize(@NotNull DataOutput2 out, @NotNull StateZ data)
			throws IOException {

				out.writeInt(data.statei);

				s.serialize(out, data.zSumBounds.lower);
				s.serialize(out, data.zSumBounds.upper);
				s.serialize(out, data.zSumDropped);

				var confEncoding = confEncodings[data.statei];
				var numPos = confSpace.states.get(data.statei).confSpace.numPos();

				out.writeInt(data.bestConfs.size());
				for (var econf : data.bestConfs) {
					assert (econf.getAssignments().length == numPos);
					for (int posi=0; posi<numPos; posi++) {
						int[] conf = econf.getAssignments();
						confEncoding.write(out, writeConf(conf[posi]));
					}
					out.writeDouble(econf.getScore());
					out.writeDouble(econf.getEnergy());
				}
			}

			@Override
			public StateZ deserialize(@NotNull DataInput2 in, int available)
			throws IOException {

				int statei = in.readInt();

				var lower = s.deserialize(in, available);
				var upper = s.deserialize(in, available);
				var dropped = s.deserialize(in, available);

				var statez = new StateZ(statei, new BigDecimalBounds(lower, upper), dropped);

				var confEncoding = confEncodings[statei];
				var numPos = confSpace.states.get(statei).confSpace.numPos();

				int numConfs = in.readInt();
				for (int i=0; i<numConfs; i++) {
					int[] conf = new int[numPos];
					for (int posi=0; posi<numPos; posi++) {
						conf[posi] = readConf(confEncoding.read(in));
					}
					double score = in.readDouble();
					double energy = in.readDouble();
					statez.bestConfs.add(new ConfSearch.EnergiedConf(conf, score, energy));
				}

				return statez;
			}

			@Override
			public int compare(StateZ a, StateZ b) {
				throw new UnsupportedOperationException();
			}
		};
	}

	public static MapDBTools.SimpleSerializer<SeqInfo> mapdbSeqInfo(MultiStateConfSpace confSpace) {
		return new MapDBTools.SimpleSerializer<>() {

			private final MapDBTools.SimpleSerializer<StateZ> s = mapdbStateZ(confSpace);

			@Override
			public void serialize(@NotNull DataOutput2 out, @NotNull SeqInfo data)
			throws IOException {
				for (int statei=0; statei<confSpace.sequencedStates.size(); statei++) {
					s.serialize(out, data.statezs[statei]);
				}
			}

			@Override
			public SeqInfo deserialize(@NotNull DataInput2 in, int available)
			throws IOException {
				SeqInfo data = new SeqInfo(confSpace);
				for (int statei=0; statei<confSpace.sequencedStates.size(); statei++) {
					data.statezs[statei] = s.deserialize(in, available);
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
		int numStates = confSpace.states.size();
		out.writeInt(numStates);
		int numSequencedStates = confSpace.sequencedStates.size();
		out.writeInt(numSequencedStates);
		var stateEncoding = IntEncoding.get(confSpace.states.size() - 1);
		out.writeByte(stateEncoding.ordinal());
		int numSeqPos = confSpace.seqSpace.positions.size();
		out.writeInt(numSeqPos);
		var seqEncoding = getSeqEncoding(confSpace.seqSpace);
		out.writeByte(seqEncoding.ordinal());
		int[] numConfPoss = confSpace.states.stream()
			.mapToInt(state -> state.confSpace.numPos())
			.toArray();
		for (var i : numConfPoss) {
			out.writeInt(i);
		}
		IntEncoding[] confEncodings = confSpace.states.stream()
			.map(state -> getConfEncoding(state.confSpace))
			.toArray(IntEncoding[]::new);
		for (var confEncoding : confEncodings) {
			out.writeByte(confEncoding.ordinal());
		}

		Consumer<StateZ> writeStatez = statez -> {
			try {

				stateEncoding.write(out, statez.statei);

				out.writeObject(statez.zSumBounds.lower);
				out.writeObject(statez.zSumBounds.upper);
				out.writeObject(statez.zSumDropped);

				var confEncoding = confEncodings[statez.statei];
				int numConfPos = numConfPoss[statez.statei];

				out.writeInt(statez.bestConfs.size());
				for (var econf : statez.bestConfs) {
					for (int posi=0; posi<numConfPos; posi++) {
						int[] conf = econf.getAssignments();
						confEncoding.write(out, writeConf(conf[posi]));
					}
					out.writeDouble(econf.getScore());
					out.writeDouble(econf.getEnergy());
				}

			} catch (IOException ex) {
				throw new RuntimeException(ex);
			}
		};

		// write the sequenced sums
		out.writeInt(op.sequencedSums.length);
		for (int sumi=0; sumi<op.sequencedSums.length; sumi++) {
			var sum = op.sequencedSums[sumi];
			for (int posi=0; posi<numSeqPos; posi++) {
				seqEncoding.write(out, writeSeq(sum.seq[posi]));
			}
			for (int statei=0; statei<numSequencedStates; statei++) {
				writeStatez.accept(sum.statezs[statei]);
			}
		}

		// write the unsequenced sums
		out.writeInt(op.unsequencedSums.length);
		for (int sumi=0; sumi<op.unsequencedSums.length; sumi++) {
			var sum = op.unsequencedSums[sumi];
			writeStatez.accept(sum.statez);
		}
	}

	public static void hazelcastSeqDBDeserializeDeNovo(ObjectDataInput in, SaveOperation op)
	throws IOException {

		// read the sequence space metadata
		int numStates = in.readInt();
		int numSequencedStates = in.readInt();
		IntEncoding stateEncoding = IntEncoding.values()[in.readByte()];
		int numPos = in.readInt();
		IntEncoding seqEncoding = IntEncoding.values()[in.readByte()];
		int[] numConfPoss = new int[numStates];
		for (int statei=0; statei<numStates; statei++) {
			numConfPoss[statei] = in.readInt();
		}
		IntEncoding[] confEncodings = new IntEncoding[numStates];
		for (int statei=0; statei<numStates; statei++) {
			confEncodings[statei] = IntEncoding.values()[in.readByte()];
		}

		Supplier<StateZ> readStatez = () -> {
			try {

				int statei = stateEncoding.read(in);

				BigDecimal lower = in.readObject();
				BigDecimal upper = in.readObject();
				BigDecimal dropped = in.readObject();

				var statez = new StateZ(
					statei,
					new BigDecimalBounds(lower, upper),
					dropped
				);

				var confEncoding = confEncodings[statei];
				int numConfPos = numConfPoss[statei];

				int numConfs = in.readInt();
				for (int i=0; i<numConfs; i++) {
					int[] conf = new int[numConfPos];
					for (int posi=0; posi<numConfPos; posi++) {
						conf[posi] = readConf(confEncoding.read(in));
					}
					double score = in.readDouble();
					double energy = in.readDouble();
					statez.bestConfs.add(new ConfSearch.EnergiedConf(conf, score, energy));
				}

				return statez;

			} catch (IOException ex) {
				throw new RuntimeException(ex);
			}
		};

		// read the sequenced sums
		op.sequencedSums = new SaveOperation.SequencedSum[in.readInt()];
		for (int sumi=0; sumi<op.sequencedSums.length; sumi++) {
			int[] seq = new int[numPos];
			for (int posi=0; posi<numPos; posi++) {
				seq[posi] = readSeq(seqEncoding.read(in));
			}
			var statezs = new StateZ[numSequencedStates];
			for (int statei=0; statei<numSequencedStates; statei++) {
				statezs[statei] = readStatez.get();
			}
			op.sequencedSums[sumi] = new SaveOperation.SequencedSum(seq, statezs);
		}

		// read the unsequenced sums
		op.unsequencedSums = new SaveOperation.UnsequencedSum[in.readInt()];
		for (int sumi=0; sumi<op.unsequencedSums.length; sumi++) {
			var statez = readStatez.get();
			op.unsequencedSums[sumi] = new SaveOperation.UnsequencedSum(statez);
		}
	}
}
