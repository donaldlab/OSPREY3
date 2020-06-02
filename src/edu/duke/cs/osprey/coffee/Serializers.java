package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.confspace.ConfSpaceIteration;
import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.tools.MapDBTools;

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

	public static MapDBTools.IntArraySerializer seq(SeqSpace seqSpace) {
		return new MapDBTools.IntArraySerializer(
			seqSpace.positions.stream()
				.flatMap(pos -> pos.resTypes.stream())
				.mapToInt(resType -> resType.index)
				.max()
				.orElse(-1),
			seqSpace.positions.size()
		);
	}
}
