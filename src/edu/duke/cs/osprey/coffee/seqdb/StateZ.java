package edu.duke.cs.osprey.coffee.seqdb;

import edu.duke.cs.osprey.tools.HashCalculator;
import edu.duke.cs.osprey.tools.Log;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;

import java.math.BigDecimal;


public class StateZ {

	public BigDecimalBounds zSumBounds;
	public BigDecimal zSumDropped;

	public StateZ(BigDecimalBounds zSumBounds, BigDecimal zSumDropped) {
		this.zSumBounds = zSumBounds;
		this.zSumDropped = zSumDropped;
	}

	public static StateZ makeUnknown() {
		return new StateZ(
			new BigDecimalBounds(BigDecimal.ZERO, MathTools.BigPositiveInfinity),
			BigDecimal.ZERO
		);
	}

	public static StateZ makeZero() {
		return new StateZ(
			BigDecimalBounds.makeZero(),
			BigDecimal.ZERO
		);
	}

	@Override
	public int hashCode() {
		return HashCalculator.combineObjHashes(zSumBounds, zSumDropped);
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof StateZ && equals((StateZ)other);
	}

	public boolean equals(StateZ other) {
		return this.zSumBounds.equals(other.zSumBounds)
			&& this.zSumDropped.equals(other.zSumDropped);
	}

	@Override
	public String toString() {
		return String.format("%s (%s)", Log.formatBigEngineering(zSumBounds), Log.formatBigEngineering(zSumDropped));
	}
}
