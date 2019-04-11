/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

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
