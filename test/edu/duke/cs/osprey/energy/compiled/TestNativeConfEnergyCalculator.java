package edu.duke.cs.osprey.energy.compiled;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.TestConfSpace;
import org.junit.Test;


public class TestNativeConfEnergyCalculator {

	@Test
	public void assign() {

		ConfSpace confSpace = TestConfSpace.Design2RL0Interface7Mut.makeCompiled().complex;
		var confEcalc = new NativeConfEnergyCalculator(confSpace, NativeConfEnergyCalculator.Precision.Double);

		assertThat(confSpace.positions.length, is(7));

		// compare the coords for a few different conformations
		int[][] confs = {
			{ 0, 0, 0, 0, 0, 0, 0 },
			{ 1, 0, 0, 0, 0, 0, 0 },
			{ 0, 1, 0, 0, 0, 0, 0 },
			{ 0, 0, 1, 0, 0, 0, 0 },
			{ 0, 0, 0, 1, 0, 0, 0 },
			{ 0, 0, 0, 0, 1, 0, 0 },
			{ 0, 0, 0, 0, 0, 1, 0 },
			{ 0, 0, 0, 0, 0, 0, 1 }
		};
		for (int[] conf : confs) {
			assertThat(confEcalc.assign(conf), is(new AssignedCoords(confEcalc.confSpace, conf).coords));
		}
	}
}
