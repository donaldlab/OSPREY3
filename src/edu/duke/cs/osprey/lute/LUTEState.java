package edu.duke.cs.osprey.lute;


import edu.duke.cs.osprey.confspace.RCTuple;

/**
 * a persistable form of all the state needed to compute conformation energies
 */
public class LUTEState {

	public final RCTuple[] tuples;
	public final double[] tupleEnergies;
	public double tupleEnergyOffset;

	public LUTEState(int numTuples) {
		tuples = new RCTuple[numTuples];
		tupleEnergies = new double[numTuples];
		tupleEnergyOffset = 0.0;
	}

	public LUTEState(LUTE.LinearSystem system) {
		this(system.tuples.size());
		set(system);
	}

	public void set(LUTE.LinearSystem system) {

		tupleEnergyOffset = system.tupleEnergyOffset;

		// copy the state from a trained LUTE
		for (int i=0; i<tuples.length; i++) {
			tuples[i] = system.tuples.get(i);
			tupleEnergies[i] = system.tupleEnergies[i];
		}
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof LUTEState && equals((LUTEState)other);
	}

	public boolean equals(LUTEState other) {

		if (this.tuples.length != other.tuples.length) {
			return false;
		}

		for (int i=0; i<tuples.length; i++) {

			if (!this.tuples[i].equals(other.tuples[i])) {
				return false;
			}

			if (this.tupleEnergies[i] != other.tupleEnergies[i]) {
				return false;
			}
		}

		return this.tupleEnergyOffset == other.tupleEnergyOffset;
	}
}
