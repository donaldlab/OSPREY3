package edu.duke.cs.osprey.confspace.compiled;


/**
 * An energetic interaction between design positions in a conformation space.
 */
public class PosInter {

	public static final int StaticPos = -1;

	/**
	 * Index in conf space of position 1, or `StaticPos` for the static atoms.
	 */
	public int posi1;

	/**
	 * Index in conf space of position 2, or `StaticPos` for the static atoms.
	 * Can be same as position 1, to signal a single (rather than pair) energy.
	 */
	public int posi2;

	public double weight;

	public PosInter(int posi1, int posi2, double weight) {
		this.posi1 = posi1;
		this.posi2 = posi2;
		this.weight = weight;
	}

	@Override
	public String toString() {
		return String.format("[%d,%d,%f]", posi1, posi2, weight);
	}
}
