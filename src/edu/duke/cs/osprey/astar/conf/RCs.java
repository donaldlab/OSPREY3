package edu.duke.cs.osprey.astar.conf;

import java.util.ArrayList;

import edu.duke.cs.osprey.pruning.PruningMatrix;

public class RCs {
	
	private int[][] unprunedRCsAtPos;
	
	public RCs(PruningMatrix pruneMat) {
		
		// pack unpruned rotamers into an efficient lookup structure
		int numPos = pruneMat.getNumPos();
        unprunedRCsAtPos = new int[numPos][];
        for (int pos=0; pos<numPos; pos++) {
        	ArrayList<Integer> srcRCs = pruneMat.unprunedRCsAtPos(pos);
        	int[] destRCs = new int[srcRCs.size()];
        	for (int i=0; i<srcRCs.size(); i++) {
        		destRCs[i] = srcRCs.get(i);
        	}
        	unprunedRCsAtPos[pos] = destRCs;
        }
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
