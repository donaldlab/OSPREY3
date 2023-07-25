package edu.duke.cs.osprey.coffee;

import static org.hamcrest.Matchers.*;
import static org.hamcrest.MatcherAssert.*;
import static edu.duke.cs.osprey.tools.Log.log;

import com.hazelcast.internal.nio.BufferObjectDataInput;
import com.hazelcast.internal.nio.BufferObjectDataOutput;
import com.hazelcast.internal.serialization.Data;
import com.hazelcast.internal.serialization.InputOutputFactory;
import com.hazelcast.internal.serialization.InternalSerializationService;
import com.hazelcast.internal.serialization.impl.ObjectDataInputStream;
import com.hazelcast.internal.serialization.impl.ObjectDataOutputStream;
import com.hazelcast.internal.serialization.impl.SerializationServiceV1;
import edu.duke.cs.osprey.coffee.nodedb.NodeIndex;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.FileTools;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;
import java.util.List;


public class DebugSerializers {

	public static void main(String[] args)
	throws Exception {

		log("reading conf space ...");
		var file = new File("path/to/confspace.ccsx");
		var confSpace = ConfSpace.fromBytes(FileTools.readFileBytes(file));
		log("done!");

		// make a multi-state conf space
		var msConfSpace = new MultiStateConfSpace.Builder("debug", confSpace)
			.build();

		// make a node out of the conformation
		var nodes = List.of(

			// wild-type conf
			new NodeIndex.Node(
				0,
				Arrays.stream(confSpace.positions)
					.mapToInt(pos ->
						Arrays.stream(pos.confs)
							.filter(c -> c.id.startsWith("wt-"))
							.mapToInt(c -> c.index)
							.findFirst()
							.orElse(0) // otherwise, pick an arbitrary conformation
					)
					.toArray(),
				new BigExp(3.2),
				new BigExp(5.4)
			),

			// all zeros
			new NodeIndex.Node(
				0,
				Arrays.stream(confSpace.positions)
					.mapToInt(pos -> 0)
					.toArray(),
				new BigExp(8.4),
				new BigExp(3.8)
			),

			// all -1
			new NodeIndex.Node(
				0,
				Arrays.stream(confSpace.positions)
					.mapToInt(pos -> -1)
					.toArray(),
				new BigExp(5.9),
				new BigExp(6.1)
			)
		);

		testIndexNode(msConfSpace, nodes);
		testHazelcastNode(msConfSpace, nodes);
	}

	private static void testIndexNode(MultiStateConfSpace confSpace, List<NodeIndex.Node> nodes) {

		var serializer = Serializers.indexNode(confSpace.states.get(0));

		var buf = ByteBuffer.allocate(serializer.bytes()*nodes.size());

		for (var node : nodes) {
			serializer.serialize(buf, node);
		}

		buf.rewind();

		for (var expNode : nodes) {
			var node = serializer.deserialize(buf);
			log("node: %s", node);
			assertThat(node, is(expNode));
		}
	}

	private static void testHazelcastNode(MultiStateConfSpace confSpace, List<NodeIndex.Node> nodes)
	throws Exception {

		var serializer = Serializers.hazelcastNode(confSpace);

		var service = SerializationServiceV1.builder()
			.withInputOutputFactory(new InputOutputFactory() {
				@Override
				public BufferObjectDataInput createInput(Data data, InternalSerializationService service) {
					throw new Error("nope");
				}

				@Override
				public BufferObjectDataInput createInput(byte[] buffer, InternalSerializationService service) {
					throw new Error("nope");
				}

				@Override
				public BufferObjectDataOutput createOutput(int size, InternalSerializationService service) {
					throw new Error("nope");
				}

				@Override
				public ByteOrder getByteOrder() {
					return ByteOrder.BIG_ENDIAN;
				}
			})
			.build();


		var baOut = new ByteArrayOutputStream();
		var out = new ObjectDataOutputStream(baOut, service);

		for (var node : nodes) {
			log("node: %s", node);
			serializer.write(out, node);
		}

		var baIn = new ByteArrayInputStream(baOut.toByteArray());
		var in = new ObjectDataInputStream(baIn, service);
		for (var expNode : nodes) {
			var node = serializer.read(in);
			log("node: %s", node);
			assertThat(node, is(expNode));
		}
	}
}
