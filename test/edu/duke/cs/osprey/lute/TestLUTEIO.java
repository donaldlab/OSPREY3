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

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.confspace.RCTuple;
import org.junit.Test;

public class TestLUTEIO {

	@Test
	public void magic() {

		byte[] bytes = LUTEIO.write(new LUTEState(0));

		assertThat((char)bytes[0], is('L'));
		assertThat((char)bytes[1], is('U'));
		assertThat((char)bytes[2], is('T'));
		assertThat((char)bytes[3], is('E'));
	}

	@Test
	public void empty() {

		LUTEState state = new LUTEState(0);

		assertThat(LUTEIO.read(LUTEIO.write(state)), is(state));
	}

	@Test
	public void populated() {

		LUTEState state = new LUTEState(3);
		state.tuples[0] = new RCTuple(4, 2);
		state.tuples[1] = new RCTuple(17, 5);
		state.tuples[2] = new RCTuple(4, 2, 17, 5);
		state.tupleEnergies[0] = 4.2;
		state.tupleEnergies[1] = 17.5;
		state.tupleEnergies[2] = -9.9;
		state.tupleEnergyOffset = 1024.8;

		assertThat(LUTEIO.read(LUTEIO.write(state)), is(state));
	}
}
