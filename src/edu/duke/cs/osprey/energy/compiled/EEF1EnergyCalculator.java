package edu.duke.cs.osprey.energy.compiled;

import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import org.joml.Vector3dc;


/**
 * Calculate energies from a compiled ConfSpace.AssignedCoords according to the EEF1 forcefield.
 *
 * WARNING: this class has internal writeable state and is NOT thread-safe!
 */
public class EEF1EnergyCalculator implements ConfSpace.EnergyCalculator {

	public static final String id = "eef1";

	private final int ffi;

	public EEF1EnergyCalculator(ConfSpace confSpace) {

		// find our forcefield index (ffi) or die trying
		ffi = confSpace.getForcefieldIndexOrThrow(id);
	}

	@Override
	public int ffi() {
		return ffi;
	}

	// TODO: this doesn't depend on the position, can we optimize it out?
	@Override
	public double calcEnergy(Vector3dc pos, double[] params) {

		double dGref = params[0];

		// easy peasy, just return the reference energy
		return dGref;
	}

	@Override
	public double calcEnergy(Vector3dc pos1, Vector3dc pos2, double[] params) {

		double vdwRadius1 = params[0];
		double lambda1 = params[1];
		double vdwRadius2 = params[2];
		double lambda2 = params[3];
		double alpha1 = params[4];
		double alpha2 = params[5];

		double r2 = pos1.distanceSquared(pos2);
		double r = Math.sqrt(r2);

		double Xij = (r - vdwRadius1)/lambda1;
		double Xji = (r - vdwRadius2)/lambda2;
		return -(alpha1*Math.exp(-Xij*Xij) + alpha2*Math.exp(-Xji*Xji))/r2;
	}
}
