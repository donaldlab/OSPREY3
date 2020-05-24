package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.confspace.ConfSpaceIteration;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.IntEncoding;
import edu.duke.cs.osprey.tools.MapDBTools;
import org.jetbrains.annotations.NotNull;
import org.mapdb.DataInput2;
import org.mapdb.DataOutput2;

import java.io.IOException;
import java.util.stream.IntStream;


public class Serializers {

	public static final MapDBTools.BigExpSerializer BigExp = new MapDBTools.BigExpSerializer();

	public static MapDBTools.IntArraySerializer conf(ConfSpaceIteration confSpace) {
		return new MapDBTools.IntArraySerializer(
			IntStream.range(0, confSpace.numPos())
				.map(posi -> confSpace.numConf(posi) - 1)
				.max()
				.orElse(-1),
			confSpace.numPos()
		);
	}

	public static class StateConfSerializer extends MapDBTools.SimpleSerializer<NodeIndex.StateConf> {

		public final MultiStateConfSpace confSpace;

		private final IntEncoding stateiEncoding;
		private final MapDBTools.IntArraySerializer[] arraySerializers;

		public StateConfSerializer(MultiStateConfSpace confSpace) {

			this.confSpace = confSpace;

			stateiEncoding = IntEncoding.get(confSpace.states.size() - 1);
			arraySerializers = confSpace.states.stream()
				.map(state -> conf(state.confSpace))
				.toArray(MapDBTools.IntArraySerializer[]::new);
		}

		@Override
		public void serialize(@NotNull DataOutput2 out, @NotNull NodeIndex.StateConf data)
		throws IOException {
			stateiEncoding.write(out, data.statei);
			arraySerializers[data.statei].serialize(out, data.conf);
		}

		@Override
		public NodeIndex.StateConf deserialize(@NotNull DataInput2 in, int available)
		throws IOException {
			int statei = stateiEncoding.read(in);
			int[] conf = arraySerializers[statei].deserialize(in, available);
			return new NodeIndex.StateConf(statei,conf);
		}
	}
	public static StateConfSerializer stateConf(MultiStateConfSpace confSpace) {
		return new StateConfSerializer(confSpace);
	}
}
