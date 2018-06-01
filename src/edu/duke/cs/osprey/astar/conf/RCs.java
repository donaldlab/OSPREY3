/*
 ** This file is part of OSPREY 3.0
 **
 ** OSPREY Protein Redesign Software Version 3.0
 ** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
 **
 ** OSPREY is free software: you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License version 2
 ** as published by the Free Software Foundation.
 **
 ** You should have received a copy of the GNU General Public License
 ** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
 **
 ** OSPREY relies on grants for its development, and since visibility
 ** in the scientific literature is essential for our success, we
 ** ask that users of OSPREY cite our papers. See the CITING_OSPREY
 ** document in this distribution for more information.
 **
 ** Contact Info:
 **    Bruce Donald
 **    Duke University
 **    Department of Computer Science
 **    Levine Science Research Center (LSRC)
 **    Durham
 **    NC 27708-0129
 **    USA
 **    e-mail: www.cs.duke.edu/brd/
 **
 ** <signature of Bruce Donald>, Mar 1, 2018
 ** Bruce Donald, Professor of Computer Science
 */

package edu.duke.cs.osprey.astar.conf;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiPredicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class RCs {

	private PruningMatrix pruneMat = null;
	private int[][] unprunedRCsAtPos;

	public RCs(SimpleConfSpace confSpace) {
		this(confSpace, (pos, resConf) -> true);
	}

	public RCs(SimpleConfSpace confSpace, BiPredicate<SimpleConfSpace.Position,SimpleConfSpace.ResidueConf> filter) {
		unprunedRCsAtPos = new int[confSpace.positions.size()][];
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			unprunedRCsAtPos[pos.index] = pos.resConfs.stream()
					.filter((resConf) -> filter.test(pos, resConf))
					.mapToInt((resConf) -> resConf.index)
					.toArray();
		}
	}

	public RCs(List<List<Integer>> rcsAtPos) {
		int n = rcsAtPos.size();
		unprunedRCsAtPos = new int[n][];
		for (int pos=0; pos<n; pos++) {
			unprunedRCsAtPos[pos] = rcsAtPos.get(pos).stream()
					.mapToInt((i) -> i)
					.toArray();
		}
	}

	public RCs(PruningMatrix pruneMat) {

		this.pruneMat = pruneMat;

		int n = pruneMat.getNumPos();

		// pack unpruned rotamers into an efficient lookup structure
		unprunedRCsAtPos = new int[n][];
		for (int pos=0; pos<n; pos++) {
			unprunedRCsAtPos[pos] = pruneMat.unprunedRCsAtPos(pos).stream()
					.mapToInt((i) -> i)
					.toArray();
		}
	}

	public RCs(RCs other, PruningMatrix pmat) {
		this(other, (pos, rc) -> !pmat.isSinglePruned(pos, rc));
		this.pruneMat = pmat;
	}

	public RCs(RCs other) {
		this(other, (pos, resConf) -> true);
	}

	public RCs(RCs other, BiPredicate<Integer,Integer> filter) {
		this.pruneMat = other.pruneMat;
		this.unprunedRCsAtPos = new int[other.unprunedRCsAtPos.length][];
		for (int i=0; i<this.unprunedRCsAtPos.length; i++) {
			final int fi = i;
			this.unprunedRCsAtPos[i] = IntStream.of(other.unprunedRCsAtPos[i])
					.filter((rc) -> filter.test(fi, rc))
					.toArray();
		}
	}

	public PruningMatrix getPruneMat() {
		return pruneMat;
	}

	public boolean hasConfs() {
		for (int[] rcs : unprunedRCsAtPos) {
			if (rcs.length > 0) {
				return true;
			}
		}
		return false;
	}

	public int getNumPos() {
		return unprunedRCsAtPos.length;
	}

	public int getNumTrivialPos() {
		int count = 0;
		for (int[] rcs : unprunedRCsAtPos) {
			if (rcs.length == 1) {
				count++;
			}
		}
		return count;
	}

	/**
	 * Counts the number of full conformations in the conformation space.
	 * The count ignores conformations disallowed by pruned RC singles,
	 * but conformations disallowed by pruned RC pairs are still counted.
	 */
	public BigInteger getNumConformations() {

		if (!hasConfs()) {
			return BigInteger.ZERO;
		}

		BigInteger count = BigInteger.ONE;
		for (int[] rcs : unprunedRCsAtPos) {
			count = count.multiply(BigInteger.valueOf(rcs.length));
		}
		return count;

		// NOTE: it's super non-trivial to account for pruned n-tuples for n > 1!
		// so we're ingoring pair pruning here entirely

		// accounting for pruned pairs might actually require doing an exponential amount of work,
		// but I'm not entirely sure about that
	}

	public int[] get(int pos) {
		return unprunedRCsAtPos[pos];
	}
	public void set(int pos, int[] rcs) {
		unprunedRCsAtPos[pos] = rcs;
	}

	public int getNum(int pos) {
		return unprunedRCsAtPos[pos].length;
	}

	public int get(int pos, int rci) {
		return unprunedRCsAtPos[pos][rci];
	}

	@Override
	public String toString() {
		return "[" + String.join(",", Arrays.stream(unprunedRCsAtPos)
				.map((int[] rcs) -> Integer.toString(rcs.length))
				.collect(Collectors.toList())
		) + "]";
	}
}