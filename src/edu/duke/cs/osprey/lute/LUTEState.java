/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

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
