package edu.duke.cs.osprey.coffee;

import com.hazelcast.nio.ObjectDataInput;
import com.hazelcast.nio.ObjectDataOutput;
import com.hazelcast.nio.serialization.StreamSerializer;
import edu.duke.cs.osprey.coffee.db.FixedIndex;
import edu.duke.cs.osprey.coffee.db.NodeIndex;
import edu.duke.cs.osprey.confspace.ConfSpaceIteration;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.IntEncoding;
import edu.duke.cs.osprey.tools.MapDBTools;

import java.io.IOException;
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
}
