package edu.duke.cs.osprey.coffee.seqdb;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.tools.Log;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;

import java.math.BigDecimal;
import java.util.*;


public class StateZ {

	public final int statei;

	public BigDecimalBounds zSumBounds;
	public BigDecimal zSumDropped;
	public final TreeSet<ConfSearch.EnergiedConf> bestConfs;

	public StateZ(int statei, BigDecimalBounds zSumBounds, BigDecimal zSumDropped) {

		this.statei = statei;

		this.zSumBounds = zSumBounds;
		this.zSumDropped = zSumDropped;

		bestConfs = new TreeSet<>(Comparator.comparing(econf -> econf.getEnergy()));
	}

	public static StateZ makeUnknown(int statei) {
		return new StateZ(
			statei,
			new BigDecimalBounds(BigDecimal.ZERO, MathTools.BigPositiveInfinity),
			BigDecimal.ZERO
		);
	}

	public static StateZ makeZero(int statei) {
		return new StateZ(
			statei,
			BigDecimalBounds.makeZero(),
			BigDecimal.ZERO
		);
	}

	public void keepBestConfs(ConfSearch.EnergiedConf econf, int num) {

		// NOTE: confs are sorted by energy, so last is the highest energy

		// if we're not tracking best confs, ignore it
		if (num <= 0) {
			return;
		}

		// if we're full and the new conf isn't better than the worst conf we have already, ignore it
		if (bestConfs.size() == num && bestConfs.comparator().compare(econf, bestConfs.last()) >= 0) {
			return;
		}

		// otherwise, add it
		bestConfs.add(econf);

		// but don't exceed the limit
		while (bestConfs.size() > num) {
			bestConfs.pollLast();
		}
	}

	@Override
	public String toString() {
		return String.format("%s (%s)", Log.formatBigEngineering(zSumBounds), Log.formatBigEngineering(zSumDropped));
	}
}
