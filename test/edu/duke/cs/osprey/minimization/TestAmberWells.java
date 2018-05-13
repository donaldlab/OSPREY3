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

package edu.duke.cs.osprey.minimization;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.Test;

import java.util.function.Consumer;

public class TestAmberWells {

	public static final double InfiniteWellEnergy = -10000.0;

	public double calcEnergy(SimpleConfSpace confSpace, RCTuple frag, Consumer<EnergyCalculator.Builder> configer) {

		EnergyCalculator.Builder builder = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setType(EnergyCalculator.Type.Cpu);
		configer.accept(builder);
		try (EnergyCalculator ecalc = builder.build()) {
			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();
			return confEcalc.calcPairEnergy(frag).energy;
		}
	}

	public void checkInfiniteWellCorrection(Strand strand, RCTuple frag) {

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		double original = calcEnergy(confSpace, frag, (b) -> {});
		assertThat(original, lessThanOrEqualTo(InfiniteWellEnergy));

		double corrected = calcEnergy(confSpace, frag, (b) -> b.setInfiniteWellEnergy(InfiniteWellEnergy));
		assertThat(corrected, greaterThan(InfiniteWellEnergy));
	}

	// these are two known cases where protein minimization falls into one of AMBER's infinitely deep energy wells

	@Test
	public void ser210vsArg193() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1gua_adj.min.pdb")).build();
		strand.flexibility.get("193").setLibraryRotamers("ARG").addWildTypeRotamers().setContinuous();
		strand.flexibility.get("210").setLibraryRotamers("SER").addWildTypeRotamers().setContinuous();
		RCTuple frag = new RCTuple(0, 28, 1, 11);

		checkInfiniteWellCorrection(strand, frag);
	}

	@Test
	public void thr210vsArg193() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1gua_adj.min.pdb")).build();
		strand.flexibility.get("193").setLibraryRotamers("ARG").addWildTypeRotamers().setContinuous();
		strand.flexibility.get("210").setLibraryRotamers("THR").addWildTypeRotamers().setContinuous();
		RCTuple frag = new RCTuple(1, 11, 0, 28);

		checkInfiniteWellCorrection(strand, frag);
	}
}
