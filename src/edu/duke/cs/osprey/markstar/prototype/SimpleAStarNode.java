package edu.duke.cs.osprey.markstar.prototype;

import java.util.Collection;

public class SimpleAStarNode {
	
	boolean isLeaf;
	public double upperBound;
	public double lowerBound;
	
	public boolean isLeaf()
	{
		return isLeaf;
	}

	public Collection<SimpleAStarNode> children() {
		return null;
	}

}
