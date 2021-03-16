package edu.duke.cs.osprey.confspace.compiled;


import java.util.Set;

/**
 * An energetic interaction between design positions in a conformation space.
 */
public class PosInter {

	public static final int StaticPos = -1;

	/**
	 * Index in conf space of position 1, or `StaticPos` for the static atoms.
	 */
	public final int posi1;

	/**
	 * Index in conf space of position 2, or `StaticPos` for the static atoms.
	 * Can be same as position 1, to signal a single (rather than pair) energy.
	 */
	public final int posi2;

	/** weight*(energy + offset) */
	public final double weight;

	/** weight*(energy + offset) */
	public final double offset;

	public PosInter(int posi1, int posi2, double weight, double offset) {
		this.posi1 = posi1;
		this.posi2 = posi2;
		this.weight = weight;
		this.offset = offset;
	}

	public boolean isIncludedIn(Set<Integer> posIndices) {
		return posIndices.contains(posi1) || posIndices.contains(posi2);
	}

	@Override
	public String toString() {
		return String.format("[%d,%d,%f,%f]", posi1, posi2, weight, offset);
	}
}
