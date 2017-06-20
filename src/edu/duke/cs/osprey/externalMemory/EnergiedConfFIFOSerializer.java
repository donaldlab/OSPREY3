package edu.duke.cs.osprey.externalMemory;

import java.nio.ByteBuffer;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.tpie.serialization.SerializingFIFOQueue;

public class EnergiedConfFIFOSerializer extends AssignmentsSerializer implements SerializingFIFOQueue.Serializer<EnergiedConf> {
	
	public EnergiedConfFIFOSerializer(RCs rcs) {
		super(rcs, Double.BYTES*2);
	}
	
	@Override
	public void serialize(EnergiedConf conf, ByteBuffer buf) {
		writeAssignments(conf.getAssignments(), buf);
		buf.putDouble(conf.getScore());
		buf.putDouble(conf.getEnergy());
	}
	
	@Override
	public EnergiedConf deserialize(ByteBuffer buf) {
		return new EnergiedConf(readAssignments(buf), buf.getDouble(), buf.getDouble());
	}
}
