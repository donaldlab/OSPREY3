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
