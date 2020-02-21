package edu.duke.cs.osprey.energy.compiled;

import java.io.DataInput;
import java.io.IOException;


/**
 * Calculate energies from a compiled ConfSpace.AssignedCoords according to the Amber forcefield.
 */
public class AmberEnergyCalculator implements EnergyCalculator {

	public static final String implementation = "amber";

	public final String id;
	public final int ffi;

	private boolean hasSettings = false;
	private boolean distanceDependentDielectric;

	public AmberEnergyCalculator(String id, int ffi) {
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
	public void readSettings(DataInput in)
	throws IOException {
		distanceDependentDielectric = in.readBoolean();
		hasSettings = true;
	}

	private void checkSettings() {
		// make sure we've read settings, or throw
		if (!hasSettings) {
			throw new IllegalStateException("haven't read settings yet");
		}
	}

	@Override
	public double calcEnergy(double r, double r2, double[] params) {

		// just in case ...
		checkSettings();

		double esQ = params[0];
		double vdwA = params[1];
		double vdwB = params[2];

		// calculate the electrostatics energy
		double es;
		if (distanceDependentDielectric) {
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
