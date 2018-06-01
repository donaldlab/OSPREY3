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
