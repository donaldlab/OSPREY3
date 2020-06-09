package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.TestConfSpace;
import edu.duke.cs.osprey.tools.FileTools;


public class TestCoffee {

	/**
	 * A 7-residue design with mutable residues on both sides
	 */
	public static MultiStateConfSpace affinity_2RL0_7mut() {
		var confSpaces = TestConfSpace.Design2RL0Interface7Mut.makeCompiled();
		return new MultiStateConfSpace
			.Builder("complex", confSpaces.complex)
			.addMutableState("A", confSpaces.chainA)
			.addMutableState("B", confSpaces.chainB)
			.build();
	}

	/**
	 * A tiny peptide inhibitor affinity design
	 */
	public static MultiStateConfSpace affinity_6ov7_1mut2flex() {
		var confSpaces = new TestConfSpace.AffinityCompiled(
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.tiny.complex.ccsx")),
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.tiny.design.ccsx")),
			ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.tiny.target.ccsx"))
		);
		return new MultiStateConfSpace
			.Builder("complex", confSpaces.complex)
			.addMutableState("design", confSpaces.chainA)
			.addUnmutableState("target", confSpaces.chainB)
			.build();
	}
}
