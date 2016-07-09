package edu.duke.cs.osprey.astar.conf;

import java.util.List;

import edu.duke.cs.osprey.pruning.PruningMatrix;

public class RCs {
	
	private PruningMatrix pruneMat;
	private int[][] unprunedRCsAtPos;
	
	public RCs(List<List<Integer>> rcsAtPos) {
		
		this.pruneMat = null;
		
		// pack the rcs into an efficient lookup structure
		int n = rcsAtPos.size();
        unprunedRCsAtPos = new int[n][];
        for (int pos=0; pos<n; pos++) {
        	List<Integer> srcRCs = rcsAtPos.get(pos);
        	int[] destRCs = new int[srcRCs.size()];
        	for (int i=0; i<srcRCs.size(); i++) {
        		destRCs[i] = srcRCs.get(i);
        	}
        	unprunedRCsAtPos[pos] = destRCs;
        }
	}
	
	public RCs(PruningMatrix pruneMat) {
		
		this.pruneMat = pruneMat;
		
		// pack unpruned rotamers into an efficient lookup structure
		int n = pruneMat.getNumPos();
        unprunedRCsAtPos = new int[n][];
        for (int pos=0; pos<n; pos++) {
        	List<Integer> srcRCs = pruneMat.unprunedRCsAtPos(pos);
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
