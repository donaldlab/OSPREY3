package edu.duke.cs.osprey.markstar.prototype;

import java.util.PriorityQueue;

public class SimpleAStarSearch {
	
	public PriorityQueue<SimpleAStarNode> heap;
	
	public SimpleConf nextBestConformation;
	private SimpleAStarNode root;
	
	public SimpleAStarSearch()
	{
		root = new SimpleAStarNode(0);
	}
	
	public SimpleAStarNode getRoot() {
		// TODO Auto-generated method stub
		return root;
	}

}
