package edu.duke.cs.osprey.kstar;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.PriorityQueue;

public class KAStarTree {

	private PriorityQueue<KAStarNode> pq = null;
	
	public KAStarTree( KSAbstract ksObj, HashMap<Integer, KSAllowedSeqs> strand2AllowedSeqs, KSCalc wt ) {
		
		// initialize KUStarNode static methods
		KAStarNode.init(ksObj, strand2AllowedSeqs, wt);
		
		pq = new PriorityQueue<KAStarNode>(strand2AllowedSeqs.get(2).getNumSeqs()*2, 
				KAStarNode.KUStarNodeComparator);
	}
	
	
	public KAStarNode poll() {
		return pq.poll();
	}
	
	
	public KAStarNode peek() {
		return pq.peek();
	}
	
	
	public int size() {
		return pq.size();
	}
	
	
	public void add( KAStarNode node ) {
		pq.add(node);
	}
	
	
	public void add( ArrayList<KAStarNode> nodes ) {
		for(KAStarNode node : nodes) add(node);
	}
}
