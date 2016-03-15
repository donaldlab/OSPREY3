package edu.duke.cs.osprey.kstar;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.PriorityQueue;

public class KUStarTree {

	private PriorityQueue<KUStarNode> pq = null;
	
	public KUStarTree( KSAbstract ksObj, HashMap<Integer, AllowedSeqs> strand2AllowedSeqs, KSCalc wt ) {
		
		// initialize KUStarNode static methods
		KUStarNode.init(ksObj, strand2AllowedSeqs, wt);
		
		pq = new PriorityQueue<KUStarNode>(strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs()/4, 
				KUStarNode.KUStarNodeComparator);
	}
	
	
	public KUStarNode poll() {
		return pq.poll();
	}
	
	
	public int size() {
		return pq.size();
	}
	
	
	public void add( KUStarNode node ) {
		pq.add(node);
	}
	
	
	public void add( ArrayList<KUStarNode> nodes ) {
		for(KUStarNode node : nodes) add(node);
	}
}
