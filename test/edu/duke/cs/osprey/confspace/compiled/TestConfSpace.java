package edu.duke.cs.osprey.confspace.compiled;


import edu.duke.cs.osprey.energy.compiled.EEF1EnergyCalculator;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;

import static edu.duke.cs.osprey.tools.Log.log;

public class TestConfSpace {

	@Test
	public void test() {

		// TEMP: read the TOML file
		String toml = FileTools.readFile("../osprey-gui/test.ccs.toml");

		ConfSpace confSpace = new ConfSpace(toml);
		EEF1EnergyCalculator ecalc = new EEF1EnergyCalculator(confSpace);

		// calculate the EEF1 energy
		log("Energy: %.6f", ecalc.calcEnergy(confSpace.assign(new int[] { 0, 0, 0 })));
	}
}
