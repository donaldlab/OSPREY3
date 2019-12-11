package edu.duke.cs.osprey.kstar.compiled;


import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;

public class TestConfEnergyCalculatorAdapter {

	@Test
	public void test() { // TODO: rename me

		ConfSpace complex = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.complex.ccs.xz"));
		ConfSpace chainA = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.A.ccs.xz"));
		ConfSpace chainG = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.G.ccs.xz"));
	}
}
