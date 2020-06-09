package edu.duke.cs.osprey.coffee.seqdb;

import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.HashCalculator;
import edu.duke.cs.osprey.tools.Log;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import edu.duke.cs.osprey.tools.Streams;

import java.math.BigDecimal;
import java.util.Arrays;


public class SeqInfo {

	public final BigDecimalBounds[] zSumBounds;

	public SeqInfo(MultiStateConfSpace confSpace) {
		this.zSumBounds = new BigDecimalBounds[confSpace.sequencedStates.size()];
	}

	public void setEmpty() {
		for (int i=0; i<zSumBounds.length; i++) {
			zSumBounds[i] = BigDecimalBounds.makeZero();
		}
	}

	public void setUnknown() {
		for (int i=0; i<zSumBounds.length; i++) {
			zSumBounds[i] = new BigDecimalBounds(BigDecimal.ZERO, MathTools.BigPositiveInfinity);
		}
	}

	public boolean isEmpty() {
		for (BigDecimalBounds b : zSumBounds) {
			if (!b.isZero()) {
				return false;
			}
		}
		return true;
	}

	public BigDecimalBounds get(MultiStateConfSpace.State state) {
		return zSumBounds[state.sequencedIndex];
	}

	@Override
	public int hashCode() {
		return HashCalculator.combineObjHashes(zSumBounds);
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof SeqInfo && equals((SeqInfo)other);
	}

	public boolean equals(SeqInfo other) {
		return Arrays.equals(this.zSumBounds, other.zSumBounds);
	}

	@Override
	public String toString() {
		return Streams.joinToString(zSumBounds, ", ", b -> Log.formatBigLn(b));
	}
}
