package edu.duke.cs.osprey.astar.conf;

import java.math.BigInteger;
import java.util.List;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class RCs {
	
	private PruningMatrix pruneMat;
	private int[][] unprunedRCsAtPos;
	private boolean hasConfs;
	private BigInteger unprunedConfsFromRCs = BigInteger.ONE;
	
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
        	unprunedConfsFromRCs = unprunedConfsFromRCs.multiply(BigInteger.valueOf(srcRCs.size()));
        }
	}
	
	public void breakDownRCSpace()
	{
		String sizeBreakdownText = "Space of all possible RC combinations: ";
		for(int i = 0; i < unprunedRCsAtPos.length-1; i++)
		{
			sizeBreakdownText+=unprunedRCsAtPos[i].length+"*";
		}
		sizeBreakdownText+=unprunedRCsAtPos[unprunedRCsAtPos.length-1].length+"="+unprunedConfsFromRCs;
		System.out.println(sizeBreakdownText);
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
        	String unprunedConfsBeforeIndex = unprunedConfsFromRCs.toString();
        	unprunedConfsFromRCs = unprunedConfsFromRCs.multiply(BigInteger.valueOf(srcRCs.size()));
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
	
	public RCs returnSubspace(RCTuple initialConf)
	{
		int[][] output = unprunedRCsAtPos.clone();
		BigInteger newSize = unprunedConfsFromRCs;
		for(int index = 0; index < initialConf.size(); index++)
		{
			int position = initialConf.pos.get(index);
			newSize = newSize.divide(BigInteger.valueOf(output[position].length));
			output[position] = new int[]{initialConf.RCs.get(index)};
		}
		return new RCs(hasConfs, pruneMat, output, newSize);
	}
	
	private RCs(boolean notEmpty, PruningMatrix pMat, int[][] unprunedRCs, BigInteger size)
	{
		unprunedRCsAtPos = unprunedRCs;
		pruneMat = pMat;
		hasConfs = notEmpty;
		unprunedConfsFromRCs = size;
	}
	
	public BigInteger unprunedConfsFromRCs()
	{
		return unprunedConfsFromRCs;
	}
}
