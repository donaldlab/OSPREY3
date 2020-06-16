package edu.duke.cs.osprey.coffee.seqdb;

import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.HashCalculator;
import edu.duke.cs.osprey.tools.Streams;

import java.util.Arrays;


public class SeqInfo {

	public final StateZ[] statezs;

	public SeqInfo(MultiStateConfSpace confSpace) {
		statezs = new StateZ[confSpace.sequencedStates.size()];
	}

	public static SeqInfo makeUnknown(MultiStateConfSpace confSpace) {
		var info = new SeqInfo(confSpace);
		for (int i = 0; i<info.statezs.length; i++) {
			info.statezs[i] = StateZ.makeUnknown();
		}
		return info;
	}

	public static SeqInfo makeZero(MultiStateConfSpace confSpace) {
		var info = new SeqInfo(confSpace);
		for (int i = 0; i<info.statezs.length; i++) {
			info.statezs[i] = StateZ.makeZero();
		}
		return info;
	}

	public StateZ get(MultiStateConfSpace.State state) {
		return statezs[state.sequencedIndex];
	}

	@Override
	public int hashCode() {
		return HashCalculator.combineObjHashes(statezs);
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof SeqInfo && equals((SeqInfo)other);
	}

	public boolean equals(SeqInfo other) {
		return Arrays.equals(this.statezs, other.statezs);
	}

	@Override
	public String toString() {
		return Streams.joinToString(statezs, ", ", s -> s.toString());
	}
}
