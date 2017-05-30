package edu.duke.cs.osprey.externalMemory;

import java.nio.ByteBuffer;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.tpie.serialization.SerializingDoublePriorityQueue;

public class EnergiedConfPrioritySerializer extends AssignmentsSerializer implements SerializingDoublePriorityQueue.Serializer<EnergiedConf> {
	
	public EnergiedConfPrioritySerializer(RCs rcs) {
		super(rcs, Double.BYTES);
	}
	
	@Override
	public double serialize(EnergiedConf conf, ByteBuffer buf) {
		writeAssignments(conf.getAssignments(), buf);
		buf.putDouble(conf.getScore());
		return conf.getEnergy();
	}
	
	@Override
	public EnergiedConf deserialize(double energy, ByteBuffer buf) {
		return new EnergiedConf(readAssignments(buf), buf.getDouble(), energy);
	}
}
