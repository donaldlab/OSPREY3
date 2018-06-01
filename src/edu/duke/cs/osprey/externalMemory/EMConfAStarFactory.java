/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.externalMemory;

import java.nio.ByteBuffer;

import edu.duke.cs.osprey.astar.conf.ConfAStarFactory;
import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.tpie.serialization.SerializingDoublePriorityQueue;

public class EMConfAStarFactory implements ConfAStarFactory {

	@Override
	public Queue<ConfAStarNode> makeQueue(RCs rcs) {
		
		Queue<EMConfAStarNode> pq = Queue.ExternalPriorityFactory.of(new NodeSerializer(rcs));
		
		// java's type system is dumb sometimes...
		Queue<? extends ConfAStarNode> q2 = (Queue<? extends ConfAStarNode>)pq;
		@SuppressWarnings("unchecked")
		Queue<ConfAStarNode> q3 = (Queue<ConfAStarNode>)q2;
		return q3;
	}

	@Override
	public ConfAStarNode makeRootNode(int numPos) {
		return new EMConfAStarNode(numPos);
	}
	
	private static class NodeSerializer extends AssignmentsSerializer implements SerializingDoublePriorityQueue.Serializer<EMConfAStarNode> {

		public NodeSerializer(RCs rcs) {
			super(rcs, Double.BYTES*2 + Integer.BYTES);
		}
		
		@Override
		public double serialize(EMConfAStarNode node, ByteBuffer buf) {
			writeAssignments(node.getConf(), buf);
			buf.putDouble(node.getGScore());
			buf.putDouble(node.getHScore());
			buf.putInt(node.getLevel());
			return node.getScore();
		}
		
		@Override
		public EMConfAStarNode deserialize(double score, ByteBuffer buf) {
			EMConfAStarNode node = new EMConfAStarNode(rcs.getNumPos());
			readAssignments(buf, node.getConf());
			node.setGScore(buf.getDouble());
			node.setHScore(buf.getDouble());
			node.setLevel(buf.getInt());
			return node;
		}
	}
}
