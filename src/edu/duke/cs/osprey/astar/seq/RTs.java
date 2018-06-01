package edu.duke.cs.osprey.astar.seq;


import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.function.BiPredicate;


/**
 * ie, residue types
 * analagous to the RCs class for residue conformations
 */
public class RTs {

	public final int numPos;

	private final int[][] typesByPos;

	/**
	 * create RTs using the full conformation space
	 */
	public RTs(SimpleConfSpace confSpace) {
		numPos = confSpace.positions.size();
		typesByPos = new int[numPos][];
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			int[] rts = new int[pos.resTypes.size()];
			for (int i=0; i<pos.resTypes.size(); i++) {
				rts[i] = i;
			}
			typesByPos[pos.index] = rts;
		}
	}

	/**
	 * select a subset of RTs using an arbitrary filter
	 */
	public RTs(RTs other, BiPredicate<Integer,Integer> filter) {
		numPos = other.typesByPos.length;
		typesByPos = new int[numPos][];
		for (int pos=0; pos<numPos; pos++) {
			final int fpos = pos;
			typesByPos[pos] = Arrays.stream(other.typesByPos[pos])
				.filter((rt) -> filter.test(fpos, rt))
				.toArray();
		}
	}

	public RTs filter(BiPredicate<Integer,Integer> filter) {
		return new RTs(this, filter);
	}

	public int numTypesAt(int pos) {
		return typesAt(pos).length;
	}

	public int numTypesAt(SimpleConfSpace.Position pos) {
		return numTypesAt(pos.index);
	}

	public int[] typesAt(int pos) {
		return typesByPos[pos];
	}

	public int[] typesAt(SimpleConfSpace.Position pos) {
		return typesAt(pos.index);
	}

	/**
	 * returns the number of full (not partial) sequences selected by this instance
	 */
	public BigInteger getNumSequences() {

		if (typesByPos.length <= 0) {
			return BigInteger.ZERO;
		}

		for (int[] rts : typesByPos) {
			if (rts.length <= 0) {
				return BigInteger.ZERO;
			}
		}

		BigInteger count = BigInteger.ONE;
		for (int[] rts : typesByPos) {
			count = count.multiply(BigInteger.valueOf(rts.length));
		}
		return count;
	}

	/**
	 * counts the number of positions with exactly one residue type
	 */
	public int getNumTrivialPos() {
		return (int)Arrays.stream(typesByPos)
			.filter(types -> types.length == 1)
			.count();
	}

	public String toString(SimpleConfSpace confSpace) {
		StringBuilder buf = new StringBuilder();
		buf.append("Residue Types:");
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			buf.append("\n\t");
			buf.append(pos.index);
			buf.append(":");
			buf.append(pos.resNum);
			buf.append("  [");
			for (int rt : typesAt(pos)) {
				buf.append(" ");
				buf.append(rt);
				buf.append(":");
				buf.append(pos.resTypes.get(rt));
			}
			buf.append(" ]");
		}
		return buf.toString();
	}
}
