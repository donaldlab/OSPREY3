package edu.duke.cs.osprey.astar.conf;

import java.util.List;

import edu.duke.cs.osprey.pruning.PruningMatrix;

public class RCs {
	
	private PruningMatrix pruneMat;
	private int[][] unprunedRCsAtPos;
	private boolean hasConfs;
	
	public RCs(List<List<Integer>> rcsAtPos) {
		
		this.pruneMat = null;
		
		int n = rcsAtPos.size();
		this.hasConfs = n > 0;
		
		// pack the rcs into an efficient lookup structure
        unprunedRCsAtPos = new int[n][];
        for (int pos=0; pos<n; pos++) {
        	List<Integer> srcRCs = rcsAtPos.get(pos);
        	
        	// no options at this pos? can't have any confs then
        	if (srcRCs.size() <= 0) {
        		hasConfs = false;
        	}
        	
        	int[] destRCs = new int[srcRCs.size()];
        	for (int i=0; i<srcRCs.size(); i++) {
        		destRCs[i] = srcRCs.get(i);
        	}
        	unprunedRCsAtPos[pos] = destRCs;
        }
	}
	
	public RCs(PruningMatrix pruneMat) {
		
		this.pruneMat = pruneMat;
		
		int n = pruneMat.getNumPos();
		this.hasConfs = n > 0;
		
		// pack unpruned rotamers into an efficient lookup structure
        unprunedRCsAtPos = new int[n][];
        for (int pos=0; pos<n; pos++) {
        	List<Integer> srcRCs = pruneMat.unprunedRCsAtPos(pos);
        	
        	// no options at this pos? can't have any confs then
        	if (srcRCs.size() <= 0) {
        		hasConfs = false;
        	}
        	
        	int[] destRCs = new int[srcRCs.size()];
        	for (int i=0; i<srcRCs.size(); i++) {
        		destRCs[i] = srcRCs.get(i);
        	}
        	unprunedRCsAtPos[pos] = destRCs;
        }
	}
	
	public PruningMatrix getPruneMat() {
		return pruneMat;
	}
	
	public boolean hasConfs() {
		return hasConfs;
	}
	
	public int getNumPos() {
		return unprunedRCsAtPos.length;
	}
	
	public int getNumTrivialPos() {
		int count = 0;
		for (int pos=0; pos<unprunedRCsAtPos.length; pos++) {
			if (unprunedRCsAtPos[pos].length == 1) {
				count++;
			}
		}
		return count;
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
}
