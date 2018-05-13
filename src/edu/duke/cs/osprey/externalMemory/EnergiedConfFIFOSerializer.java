/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

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
