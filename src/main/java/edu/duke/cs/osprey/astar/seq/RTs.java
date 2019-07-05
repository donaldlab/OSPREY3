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

package edu.duke.cs.osprey.astar.seq;


import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;
import edu.duke.cs.osprey.confspace.SeqSpace;

import java.math.BigInteger;
import java.util.Arrays;


/**
 * ie, residue types
 * analagous to the RCs class for residue conformations
 */
public class RTs {

	public final int numPos;

	private final int[] wildTypesByPos;
	private final int[][] indicesByPos;

	/**
	 * create RTs using the full sequence space
	 */
	public RTs(SeqSpace seqSpace) {
		numPos = seqSpace.positions.size();
		wildTypesByPos = new int[numPos];
		indicesByPos = new int[numPos][];
		for (SeqSpace.Position pos : seqSpace.positions) {
			wildTypesByPos[pos.index] = pos.wildType != null ? pos.wildType.index : -1;
			int[] indices = new int[pos.resTypes.size()];
			for (int i=0; i<pos.resTypes.size(); i++) {
				indices[i] = i;
			}
			indicesByPos[pos.index] = indices;
		}
	}

	public interface Filter {
		boolean test(int pos, int rt);
	}

	/**
	 * select a subset of RTs using an arbitrary filter
	 */
	public RTs(RTs other, Filter filter) {
		numPos = other.numPos;
		wildTypesByPos = other.wildTypesByPos.clone();
		indicesByPos = new int[numPos][];
		for (int mpos = 0; mpos< numPos; mpos++) {
			final int fmpos = mpos;
			indicesByPos[mpos] = Arrays.stream(other.indicesByPos[mpos])
				.filter((rt) -> filter.test(fmpos, rt))
				.toArray();
		}
	}

	public RTs filter(Filter filter) {
		return new RTs(this, filter);
	}

	public int numTypesAt(int pos) {
		return indicesAt(pos).length;
	}

	public int[] indicesAt(int pos) {
		return indicesByPos[pos];
	}

	public int wildTypeAt(int pos) {
		return wildTypesByPos[pos];
	}

	/**
	 * returns the number of full (not partial) sequences selected by this instance
	 */
	public BigInteger getNumSequences() {

		if (numPos <= 0) {
			return BigInteger.ZERO;
		}

		for (int[] rts : indicesByPos) {
			if (rts.length <= 0) {
				return BigInteger.ZERO;
			}
		}

		BigInteger count = BigInteger.ONE;
		for (int[] rts : indicesByPos) {
			count = count.multiply(BigInteger.valueOf(rts.length));
		}
		return count;
	}

	/**
	 * counts the number of positions with exactly one residue type
	 */
	public int getNumTrivialPos() {
		return (int)Arrays.stream(indicesByPos)
			.filter(indices -> indices.length == 1)
			.count();
	}

	public int getNumMutations(SeqAStarNode.Assignments assignments) {
		int count = 0;
		for (int i=0; i<assignments.numAssigned; i++) {
			int pos = assignments.assignedPos[i];
			int rt = assignments.assignedRTs[i];
			if (rt != wildTypesByPos[pos]) {
				count++;
			}
		}
		return count;
	}
}
