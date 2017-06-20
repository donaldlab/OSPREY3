package edu.duke.cs.osprey.externalMemory;

import java.util.Arrays;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;

public class EMConfAStarNode implements ConfAStarNode {
	
	public static final int NotAssigned = -1;
	
	// NOTE: we'll only keep a few nodes in memory at once,
	// so there's less pressure here to optimize space.
	// also, nodes have to be completely self-contained
	// meaning, they can't refer to other nodes for information (eg, links)
	private int[] assignments;
	private int level;
	private double gscore;
	private double hscore;
	
	public EMConfAStarNode(int numPos) {
		assignments = new int[numPos];
		Arrays.fill(assignments, NotAssigned);
		level = 0;
		gscore = 0.0;
		hscore = 0.0;
	}
	
	public EMConfAStarNode(EMConfAStarNode other) {
		assignments = Arrays.copyOf(other.assignments, other.assignments.length);
		level = other.level;
		gscore = other.gscore;
		hscore = other.hscore;
	}

	@Override
	public EMConfAStarNode assign(int pos, int rc) {
		EMConfAStarNode other = new EMConfAStarNode(this);
		other.assignments[pos] = rc;
		other.level++;
		return other;
	}

	@Override
	public double getGScore() {
		return gscore;
	}

	@Override
	public void setGScore(double val) {
		gscore = val;
	}

	@Override
	public double getHScore() {
		return hscore;
	}

	@Override
	public void setHScore(double val) {
		hscore = val;
	}

	@Override
	public int getLevel() {
		return level;
	}
	
	public void setLevel(int val) {
		level = val;
	}

	@Override
	public void getConf(int[] out) {
		System.arraycopy(assignments, 0, out, 0, assignments.length);
	}
	
	public int[] getConf() {
		return assignments;
	}
	
	@Override
	public void index(ConfIndex index) {
		
		// is this node already indexed?
		if (index.node == this) {
			return;
		}
		index.node = this;
		
		// copy values and references to the stack for speed
		int n = index.numPos;
		int numDefined = 0;
		int[] dpos = index.definedPos;
		int[] rcs = index.definedRCs;
		int numUndefined = 0;
		int[] upos = index.undefinedPos;
		
		// split conformation into defined and undefined positions
		numDefined = 0;
		numUndefined = 0;
		for (int pos=0; pos<n; pos++) {
			int rc = assignments[pos];
			
			if (rc == EMConfAStarNode.NotAssigned) {
				upos[numUndefined] = pos;
				numUndefined++;
			} else {
				dpos[numDefined] = pos;
				rcs[numDefined] = rc;
				numDefined++;
			}
		}
		
		// copy stack vars back to the index
		index.numDefined = numDefined;
		index.numUndefined = numUndefined;
	}
}
