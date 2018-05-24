/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
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
