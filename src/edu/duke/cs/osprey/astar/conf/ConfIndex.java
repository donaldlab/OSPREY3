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
	
	public boolean isDefined(int pos) {
		for (int i=0; i<numDefined; i++) {
			if (definedPos[i] == pos) {
				return true;
			}
		}
		return false;
	}
	
	public boolean isUndefined(int pos) {
		for (int i=0; i<numUndefined; i++) {
			if (undefinedPos[i] == pos) {
				return true;
			}
		}
		return false;
	}
	
	public ConfIndex assign(int nextPos, int nextRc) {
		
		ConfIndex other = new ConfIndex(numPos);
		
		// the next pos should be undefined (and not defined)
		assert (this.isUndefined(nextPos));
		assert (!this.isDefined(nextPos));
		
		// copy from the other index
		other.numDefined = this.numDefined + 1;
		other.numUndefined = this.numUndefined - 1;
		
		// update defined side
		boolean isInserted = false;
		for (int i=0; i<this.numDefined; i++) {
			int pos = this.definedPos[i];
			int rc = this.definedRCs[i];
			if (nextPos > pos) {
				other.definedPos[i] = pos;
				other.definedRCs[i] = rc;
			} else {
				
				if (!isInserted) {
					other.definedPos[i] = nextPos;
					other.definedRCs[i] = nextRc;
					isInserted = true;
				}
				
				other.definedPos[i+1] = pos;
				other.definedRCs[i+1] = rc;
			}
		}
		if (!isInserted) {
			other.definedPos[this.numDefined] = nextPos;
			other.definedRCs[this.numDefined] = nextRc;
		}
		
		// update undefined side
		int j = 0;
		for (int i=0; i<this.numUndefined; i++) {
			int pos = this.undefinedPos[i];
			if (pos != nextPos) {
				other.undefinedPos[j++] = pos;
			}
		}
		
		// init defaults for things we won't copy
		other.node = null;
		
		return other;
	}
}
