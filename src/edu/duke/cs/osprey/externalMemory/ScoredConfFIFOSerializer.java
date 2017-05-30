package edu.duke.cs.osprey.externalMemory;

import java.nio.ByteBuffer;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.tpie.serialization.SerializingFIFOQueue;

public class ScoredConfFIFOSerializer extends AssignmentsSerializer implements SerializingFIFOQueue.Serializer<ScoredConf> {
	
	public ScoredConfFIFOSerializer(RCs rcs) {
		super(rcs, Double.BYTES);
	}
	
	@Override
	public void serialize(ScoredConf conf, ByteBuffer buf) {
		writeAssignments(conf.getAssignments(), buf);
		buf.putDouble(conf.getScore());
	}
	
	@Override
	public ScoredConf deserialize(ByteBuffer buf) {
		return new ScoredConf(readAssignments(buf), buf.getDouble());
	}
}
