package edu.duke.cs.osprey.energy.compiled;

import java.io.DataInput;
import java.io.IOException;


/**
 * Calculate energies from a compiled ConfSpace.AssignedCoords according to the Amber forcefield.
 */
public class AmberEnergyCalculator implements EnergyCalculator {

	public static final Type type = Type.Amber;

	public final String id;
	public final int ffi;

	public static class Settings {
		boolean distanceDependentDielectric;
	}

	public Settings settings = null;

	public AmberEnergyCalculator(String id, int ffi) {
		this.id = id;
		this.ffi = ffi;
	}

	@Override
	public String id() {
		return id;
	}

	@Override
	public Type type() {
		return type;
	}

	@Override
	public int ffi() {
		return ffi;
	}

	@Override
	public void readSettings(DataInput in)
	throws IOException {
		settings = new Settings();
		settings.distanceDependentDielectric = in.readBoolean();
	}

	@Override
	public double calcEnergy(double r, double r2, double[] params) {

		double esQ = params[0];
		double vdwA = params[1];
		double vdwB = params[2];

		// calculate the electrostatics energy
		double es;
		if (settings.distanceDependentDielectric) {
			es = esQ/r2;
		} else {
			es = esQ/r;
		}

		// calculate the van der Waals energy
		double r6 = r2*r2*r2;
		double r12 = r6*r6;
		double vdw = vdwA/r12 - vdwB/r6;

		return es + vdw;
	}
}
