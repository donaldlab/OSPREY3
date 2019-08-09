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


import java.util.Arrays;

public class ConfIndex {
	
	public final int numPos;
	
	public ConfAStarNode node;
	public int numDefined;
	public final int[] definedPos;
	public final int[] definedRCs;
	public int numUndefined;
	public final int[] undefinedPos;
	
	public ConfIndex(int numPos) {
		this.numPos = numPos;
		this.node = null;
        this.numDefined = 0;
        this.definedPos = new int[numPos];
        this.definedRCs = new int[numPos];
        this.numUndefined = 0;
        this.undefinedPos = new int[numPos];
	}

	public ConfIndex(ConfIndex other) {
		this.numPos = other.numPos;
		this.numDefined = other.numDefined;
		this.definedPos = other.definedPos.clone();
		this.definedRCs = other.definedRCs.clone();
		this.numUndefined = other.numUndefined;
		this.undefinedPos = other.undefinedPos.clone();
		this.node = null;
	}

	public boolean isFullyDefined() {
		return numDefined == numPos;
	}
	
	public boolean isDefined(int pos) {
		return findDefined(pos) >= 0;
	}

	public int findDefined(int pos) {
		for (int i=0; i<numDefined; i++) {
			if (definedPos[i] == pos) {
				return i;
			}
		}
		return -1;
	}

	public boolean isUndefined(int pos) {
		return findUndefined(pos) >= 0;
	}

	public int findUndefined(int pos) {
		for (int i=0; i<numUndefined; i++) {
			if (undefinedPos[i] == pos) {
				return i;
			}
		}
		return -1;
	}
	
	public ConfIndex assign(int nextPos, int nextRc) {
		
		ConfIndex other = new ConfIndex(this);
		other.assignInPlace(nextPos, nextRc);
		
		// init defaults for things we won't copy
		other.node = null;
		
		return other;
	}

	public void assignInPlace(int pos, int rc) {

		// update defined side
		int insertIndex = Arrays.binarySearch(definedPos, 0, numDefined, pos);
		if (insertIndex >= 0) {
			throw new IllegalArgumentException("pos " + pos + " already assigned");
		}
		insertIndex = -insertIndex - 1;
		for (int i=numDefined; i>insertIndex; i--) {
			definedPos[i] = definedPos[i-1];
			definedRCs[i] = definedRCs[i-1];
		}
		definedPos[insertIndex] = pos;
		definedRCs[insertIndex] = rc;
		numDefined++;

		updateUndefined();
	}

	public ConfIndex unassign(int pos) {

		ConfIndex other = new ConfIndex(this);
		other.unassignInPlace(pos);

		// init defaults for things we won't copy
		other.node = null;

		return other;
	}

	public void unassignInPlace(int pos) {

		// update defined side
		int removeIndex = Arrays.binarySearch(definedPos, 0, numDefined, pos);
		if (removeIndex < 0) {
			throw new IllegalArgumentException("pos " + pos + " not assigned");
		}
		numDefined--;
		for (int i=removeIndex; i<numDefined; i++) {
			definedPos[i] = definedPos[i+1];
			definedRCs[i] = definedRCs[i+1];
		}

		updateUndefined();
	}

	/**
	 * ensures assigned and unassigned positions are sorted in increasing order
	 */
	public void sortDefined() {

		// these arrays are typically tiny (<20), so insertion sort is very efficient here
		// also, we have to sort the assigned positions and res confs together,
		// so we can't use a library sort =(

		// sort the defined side
		for (int i=1; i<numDefined; i++) {

			int tempPos = definedPos[i];
			int tempRT = definedRCs[i];

			int j;
			for (j=i; j>=1 && tempPos < definedPos[j-1]; j--) {
				definedPos[j] = definedPos[j-1];
				definedRCs[j] = definedRCs[j-1];
			}
			definedPos[j] = tempPos;
			definedRCs[j] = tempRT;
		}
	}

	/**
	 * Populates the unassigned positions, based on what's not assigned
	 * defined positions must be sorted
	 */
	public void updateUndefined() {

		numUndefined = 0;

		int i = 0;
		for (int pos=0; pos<numPos; pos++) {

			// does this pos match the next defined pos?
			if (i < numDefined && pos == definedPos[i]) {

				// yup, skip this pos
				i++;

			} else {

				// nope, it's undefined, append it
				undefinedPos[numUndefined++] = pos;
			}
		}
	}

	@Override
	public String toString() {
		StringBuilder buf = new StringBuilder();
		buf.append('[');
		for (int i=0; i<numDefined; i++) {
			if (i > 0) {
				buf.append(", ");
			}
			buf.append(definedPos[i]);
			buf.append('=');
			buf.append(definedRCs[i]);
		}
		buf.append(']');
		return buf.toString();
	}
}
