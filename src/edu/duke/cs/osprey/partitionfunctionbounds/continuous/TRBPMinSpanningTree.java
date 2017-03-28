package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.PriorityQueue;

import com.sun.javafx.geom.Edge;

public class TRBPMinSpanningTree {

	public double[][] getMinSpanningTree(double[][] edgeWeights) { 
		int numNodes = edgeWeights.length;
		UnionFind elems = new UnionFind(numNodes); 
		
		PriorityQueue<Edge> edgeQueue = new PriorityQueue<>();
		for (int i=0; i<numNodes; i++) { 
			for (int j=0; j<numNodes; j++) { 
				if (i!=j) { edgeQueue.add(new Edge(i, j, edgeWeights[i][j])); }
			}
		}
		
		ArrayList<Edge> minSpanningTree = new ArrayList<>();
		while (!edgeQueue.isEmpty() && minSpanningTree.size()<numNodes-1) {
			Edge e = edgeQueue.poll();
			int n1 = e.node1;
			int n2 = e.node2;
			if (elems.connected(n1, n2)) { continue; } 
			elems.union(n1, n2);
			minSpanningTree.add(e);
		}
		
		double[][] mstVector = new double[numNodes][numNodes];
		for (double[] r : mstVector) { for (double d : r) { d = 0; } } 
		for (Edge e : minSpanningTree) { 
			mstVector[e.node1][e.node2] = 1;
			mstVector[e.node2][e.node1] = 1;
		}
		
		return mstVector;

	}

	private class Edge implements Comparable { 
		int node1; 
		int node2;
		double weight;

		public Edge (int n1, int n2, double w) { 
			node1 = n1;
			node2 = n2;
			weight = w;
		}

		public int compareTo(Object o) {
			Edge e = (Edge) o;
			return Double.compare(this.weight, e.weight);
		}



	}
	
	private class UnionFind { 
		int[] parent;
		int[] size;
		int count;
		
		public UnionFind(int numElements) { 
			parent = new int[numElements];
			for (int i=0; i<parent.length; i++) { parent[i] = i; }
			size = new int[numElements];
			count = numElements;
		}
		
		public int countSets() { return count; } 
		
		public boolean connected(int x, int y) { 
			return parent[x] == parent[y];
		}
		
		public int find(int n) { 
			while (n != parent[n]) { n = parent[n]; }
			return n;
		}
		
		public void union(int x, int y) { 
			int rx = find(x); 
			int ry = find(y); 
			if (rx==ry) { return ; } 
			
			if (size[rx] < size[ry]) { 
				parent[rx] = ry;
				size[ry] += size[rx];
			} else { 
				parent[ry] = rx;
				size[rx] += size[ry];
			}
			count--;
		}
	}
	
}
