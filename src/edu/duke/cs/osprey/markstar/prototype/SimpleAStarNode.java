package edu.duke.cs.osprey.markstar.prototype;

import java.util.Collection;
import java.util.HashSet;

public class SimpleAStarNode {
	
	boolean isLeaf;
	public double upperBound;
	public double lowerBound;
	
	private Collection<SimpleAStarNode> children = null;
	static final int max_depth = 2;
	static final int branch_factor = 2;
	private int depth;
	
	public SimpleAStarNode(int d)
	{
		this(d, null);
	}
	
	public SimpleAStarNode(int d, SimpleConf simpleConf) {
		depth = d;
		isLeaf = depth >= max_depth;
		System.out.println("Made node with depth "+depth);
	}

	public boolean isLeaf()
	{
		return isLeaf;
	}

	public Collection<SimpleAStarNode> children() {
		if(!isLeaf)
			return generateChildren();
		else return null;
	}

	private Collection<SimpleAStarNode> generateChildren() {
		if(children != null)
			return children;
		children = new HashSet<SimpleAStarNode>();
		for(int i = 0; i < 4; i++)
		{
			children.add(new SimpleAStarNode(this.depth+1,new SimpleConf(i)));
			System.out.println("New node depth "+(depth+1));
		}
		return children;
	}
	
	public String toString()
	{
		return "";
	}
	
	private String toTreeString()
	{
		return "";
	}

}
