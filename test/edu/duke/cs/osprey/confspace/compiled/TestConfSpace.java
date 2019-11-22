package edu.duke.cs.osprey.confspace.compiled;


import edu.duke.cs.osprey.energy.compiled.EnergyCalculator;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;

import static edu.duke.cs.osprey.tools.Log.log;

public class TestConfSpace {

	@Test
	public void test() {

		// TEMP: read the TOML file
		String toml = FileTools.readFile("../osprey-gui/test.ccs.toml");

		ConfSpace confSpace = new ConfSpace(toml);
		EnergyCalculator ecalcAmber = confSpace.ecalcs[0];
		EnergyCalculator ecalcEef1 = confSpace.ecalcs[1];

		// calculate the EEF1 energy
		ConfSpace.AssignedCoords conf = confSpace.assign(new int[]{0, 0, 0});
		log("Energy Amber: %.6f", ecalcAmber.calcEnergy(conf));
		log("Energy EEF1:  %.6f", ecalcEef1.calcEnergy(conf));
	}
}
