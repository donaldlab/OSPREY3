package edu.duke.cs.osprey.astar.seq;


import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.math.BigInteger;
import java.util.Arrays;


/**
 * ie, residue types
 * analagous to the RCs class for residue conformations
 */
public class RTs {

	public final int numMutablePos;

	private final String[] resNums;
	private final String[][] namesByPos;
	private final int[][] indicesByPos;

	/**
	 * create RTs using the full conformation space
	 */
	public RTs(SimpleConfSpace confSpace) {
		numMutablePos = confSpace.mutablePositions.size();
		resNums = new String[numMutablePos];
		namesByPos = new String[numMutablePos][];
		indicesByPos = new int[numMutablePos][];
		for (SimpleConfSpace.Position pos : confSpace.mutablePositions) {
			resNums[pos.mindex] = pos.resNum;
			String[] names = new String[pos.resTypes.size()];
			int[] indices = new int[pos.resTypes.size()];
			for (int i=0; i<pos.resTypes.size(); i++) {
				indices[i] = i;
				names[i] = pos.resTypes.get(i);
			}
			namesByPos[pos.index] = names;
			indicesByPos[pos.index] = indices;
		}
	}

	public interface Filter {
		boolean test(int mposIndex, String resNum, int rtIndex, String rtName);
	}

	/**
	 * select a subset of RTs using an arbitrary filter
	 */
	public RTs(RTs other, Filter filter) {
		numMutablePos = other.numMutablePos;
		resNums = other.resNums.clone();
		namesByPos = new String[numMutablePos][];
		indicesByPos = new int[numMutablePos][];
		for (int mpos = 0; mpos<numMutablePos; mpos++) {
			final int fmpos = mpos;
			indicesByPos[mpos] = Arrays.stream(other.indicesByPos[mpos])
				.filter((rt) -> filter.test(fmpos, resNums[fmpos], rt, other.namesByPos[fmpos][rt]))
				.toArray();
			namesByPos[mpos] = Arrays.stream(indicesByPos[mpos])
				.mapToObj(rt -> other.namesByPos[fmpos][rt])
				.toArray(size -> new String[size]);
		}
	}

	public RTs filter(Filter filter) {
		return new RTs(this, filter);
	}

	public int numTypesAt(int mpos) {
		return indicesAt(mpos).length;
	}

	public int numTypesAt(SimpleConfSpace.Position mpos) {
		return numTypesAt(mpos.mutableIndexOrThrow());
	}

	public int[] indicesAt(int mpos) {
		return indicesByPos[mpos];
	}

	public int[] indicesAt(SimpleConfSpace.Position mpos) {
		return indicesAt(mpos.mutableIndexOrThrow());
	}

	public String getName(int mpos, int rt) {
		return namesByPos[mpos][rt];
	}


	/**
	 * returns the number of full (not partial) sequences selected by this instance
	 */
	public BigInteger getNumSequences() {

		if (numMutablePos <= 0) {
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

	@Override
	public String toString() {
		StringBuilder buf = new StringBuilder();
		buf.append("Residue Types:");
		for (int mpos=0; mpos<numMutablePos; mpos++) {
			buf.append("\n\t");
			buf.append(mpos);
			buf.append(":");
			buf.append(resNums[mpos]);
			buf.append("  [");
			for (int i=0; i<numTypesAt(mpos); i++) {
				buf.append(" ");
				buf.append(indicesByPos[mpos][i]);
				buf.append(":");
				buf.append(namesByPos[mpos][i]);
			}
			buf.append(" ]");
		}
		return buf.toString();
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof RTs && equals((RTs)other);
	}

	public boolean equals(RTs other) {

		if (this.numMutablePos != other.numMutablePos) {
			return false;
		}

		if (!Arrays.equals(this.resNums, other.resNums)) {
			return false;
		}

		for (int mpos=0; mpos<numMutablePos; mpos++) {
			if (!Arrays.equals(this.indicesByPos[mpos], other.indicesByPos[mpos])) {
				return false;
			}
			if (!Arrays.equals(this.namesByPos[mpos], other.namesByPos[mpos])) {
				return false;
			}
		}

		return true;
	}
}
