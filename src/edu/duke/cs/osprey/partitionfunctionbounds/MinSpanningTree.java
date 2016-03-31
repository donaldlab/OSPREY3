/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

import edu.duke.cs.osprey.tools.CreateMatrix;
import java.util.ArrayList;
import java.util.PriorityQueue;

/**
 *
 * @author hmn5
 */
public class MinSpanningTree {

    ArrayList<Edge> mst;
    public int[][] mstVector;//0-1 vector, where 1 indicates that edge 
    //between nodes i,j is in MST
    //mstVector is implemented such that each entry i, point to all nodes j s.t. j<i;
    ///this is just like the pairwise energy matrix

    public MinSpanningTree(double[][] edgeWeights, boolean[][] interactionGraph) {
        int numNodes = interactionGraph.length;
        UnionFind uf = new UnionFind(numNodes);
        PriorityQueue<Edge> pq = new PriorityQueue<>();
        //Create edges and add to queue
        for (int nodeI = 0; nodeI < numNodes; nodeI++) {
            for (int nodeJ = 0; nodeJ < nodeI; nodeJ++) {
                if (interactionGraph[nodeI][nodeJ]) {
                    Edge edge = new Edge(nodeI, nodeJ, edgeWeights[nodeI][nodeJ]);
                    pq.add(edge);
                }
            }
        }

        //begin MST algorithm
        mst = new ArrayList<>();
        while (!pq.isEmpty() && mst.size() < numNodes - 1) {
            Edge e = pq.poll();
            int v1 = e.node1;
            int v2 = e.node2;
            if (uf.connected(v1, v2)) {
                continue;
            }
            uf.union(v1, v2);
            mst.add(e);
        }

        this.mstVector = new int[numNodes][];
        //initialize to 0 for all entries
        for (int i = 0; i < numNodes; i++) {
            int[] nodesBelowI = new int[i];
            for (int j = 0; j < i; j++) {
                nodesBelowI[j] = 0;
            }
            this.mstVector[i] = nodesBelowI;
        }
        
        for (Edge e : this.mst) {
            int v1 = e.node1;
            int v2 = e.node2;
            if (v1 < v2) {
                this.mstVector[v2][v1] = 1;
            } else {
                this.mstVector[v1][v2] = 1;
            }
        }
    }

    private class Edge implements Comparable {

        int node1;
        int node2;
        double weight;

        public Edge(int node1, int node2, double weight) {
            this.node1 = node1;
            this.node2 = node2;
            this.weight = weight;
        }

        @Override
        public int compareTo(Object o) {
            Edge edge2 = (Edge) o;//we can only compare to other edges, and expect no other cases
            return Double.valueOf(this.weight).compareTo(edge2.weight);
        }
    }
}
