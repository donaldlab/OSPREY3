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
import java.util.PriorityQueue;
import java.util.Set;

public class RecursiveKStarBound {
	
	private SimpleAStarSearch baseSearchTree;
	private PriorityQueue<EpsilonNode> epsilonSortedHeap;
	private double targetEpsilon = 0;
	
	public RecursiveKStarBound(SimpleAStarSearch confAStar, double epsilon)
	{
		baseSearchTree = confAStar;
	}
	
	public void computeBound()
	{
		double currentEpsilon = boundSearchNode(baseSearchTree.getRoot());
		
		while (currentEpsilon > targetEpsilon)
		{
			EpsilonNode nextNode = epsilonSortedHeap.poll();
			if(nextNode.isLeaf())
				refineSingleSequencebound(nextNode);
			else
				branchAndBound(nextNode);
		}
	}

	private void branchAndBound(EpsilonNode nextNode) {
		for(SimpleAStarNode childNode: expand(nextNode.searchTreeNode))
			epsilonSortedHeap.add(new EpsilonNode(childNode, boundSearchNode(childNode)));
		
	}

	private Collection<SimpleAStarNode> expand(SimpleAStarNode searchTreeNode) {
		return searchTreeNode.children();
	}

	private void refineSingleSequencebound(EpsilonNode nextNode) {
		// TODO Auto-generated method stub
		
	}

	private double boundSearchNode(SimpleAStarNode node) {
		double nodeBound = node.upperBound - node.lowerBound;
		if(node.isLeaf() || nodeBound < targetEpsilon)
			return nodeBound;
		else 
		{
			double sum = 0;
			for(SimpleAStarNode childNode: expand(node))
				sum += boundSearchNode(childNode);
			return sum;
		}
				
	}
	
	private double recursivelyBoundSearchNode(SimpleAStarNode curNode)
	{
		return 0;
	}
	
	private class EpsilonNode implements Comparable<EpsilonNode> 
	{
		
		public SimpleAStarNode searchTreeNode;
		boolean isLeaf;
		private double epsilonBound = 0;
		
		public EpsilonNode(SimpleAStarNode childNode, double boundSearchNode) {
			searchTreeNode = childNode;
		}

		public boolean isLeaf()
		{
			return isLeaf;
		}

		@Override
		public int compareTo(EpsilonNode o) {
			return Double.compare(epsilonBound,o.epsilonBound);
		}
		
	}

}
