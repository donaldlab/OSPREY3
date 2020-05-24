package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.compiled.TestConfSpace;


public class TestCoffee {

	public static MultiStateConfSpace affinity_2RL0_7mut() {
		var confSpaces = TestConfSpace.Design2RL0Interface7Mut.makeCompiled();
		return new MultiStateConfSpace
			.Builder("complex", confSpaces.complex)
			.addMutableState("design", confSpaces.chainA)
			.addMutableState("target", confSpaces.chainB)
			.build();
	}
}
