package edu.duke.cs.osprey.markstar.prototype;

import java.util.Collection;
import java.util.PriorityQueue;
import java.util.Set;

public class RecursiveKStarBound {
	
	private SimpleAStarSearch baseSearchTree;
	private PriorityQueue<EpsilonNode> epsilonSortedHeap;
	private double targetEpsilon = 0;
	
	public RecursiveKStarBound(SimpleAStarSearch confAStar)
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
