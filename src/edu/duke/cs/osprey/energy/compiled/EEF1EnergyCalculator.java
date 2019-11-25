package edu.duke.cs.osprey.energy.compiled;


/**
 * Calculate energies from a compiled ConfSpace.AssignedCoords according to the EEF1 forcefield.
 */
public class EEF1EnergyCalculator implements EnergyCalculator {

	public static final String implementation = "eef1";

	// Only need to count interactions with 9 A.
	// Farther interactions are already counted in the pre-calculated internal energy.
	private static final double cutoff = 9.0;

	public final String id;
	public final int ffi;

	public EEF1EnergyCalculator(String id, int ffi) {
		this.id = id;
		this.ffi = ffi;
	}

	@Override
	public String id() {
		return id;
	}

	@Override
	public int ffi() {
		return ffi;
	}

	@Override
	public double calcEnergy(double r, double r2, double[] params) {

		double vdwRadius1 = params[0];
		double lambda1 = params[1];
		double vdwRadius2 = params[2];
		double lambda2 = params[3];
		double alpha1 = params[4];
		double alpha2 = params[5];

		if (r <= cutoff) {
			double Xij = (r - vdwRadius1)/lambda1;
			double Xji = (r - vdwRadius2)/lambda2;
			return -(alpha1*Math.exp(-Xij*Xij) + alpha2*Math.exp(-Xji*Xji))/r2;
		} else {
			return 0.0;
		}
	}
}
