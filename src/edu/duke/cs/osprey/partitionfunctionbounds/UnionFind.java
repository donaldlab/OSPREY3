/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

/**
 * This provides an implementation of UnionFind with lgN union and lgN find
 * based on the WeightedUnionFind implementation in Sedgewick and Wayne
 * 
 * @author hmn5
 */
public class UnionFind {
    
    private int[] parent;
    private int[] size; 
    private int count;
    
    public UnionFind(int numElements){
        this.parent = new int[numElements];
        this.size = new int[numElements];
        this.count = numElements;
        for (int i=0; i<numElements; i++){
            this.parent[i] = i;
            this.size[i] = 1;
        }
    }
    
    public int count(){
        return count;
    }
    
    public int find(int p){
        validate(p);
        while (p != this.parent[p]){
            p = parent[p];
        }
        return p;
    }
    
    public boolean connected(int p, int q){
        return find(p) == find(q);
    }
    
    public void union(int p, int q){
        int rootP = find(p);
        int rootQ = find(q);
        if (rootP == rootQ) return;
        
        //make smaller root point to larger root
        if (size[rootP] < size[rootQ]){
            parent[rootP] = rootQ;
            size[rootQ] += size[rootP];
        }
        else{
            parent[rootQ] = rootP;
            size[rootP] += size[rootQ];
        }
        count--;
    }
    
    private void validate(int p){
        int numElements = this.parent.length;
        if (p < 0 || p >= numElements){
            throw new IndexOutOfBoundsException("index "+ p + " is not between 0 and " + (numElements-1));
        }
    }
}
