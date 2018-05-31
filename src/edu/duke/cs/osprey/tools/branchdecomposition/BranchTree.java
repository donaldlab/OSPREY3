package edu.duke.cs.osprey.tools.branchdecomposition;
import java.io.*;
import java.util.LinkedHashSet;

public class BranchTree implements Serializable {

	private BranchNode bn[] = null; //the nodes in the tree	
	private int numNodes = 0; //the number of nodes in the tree
	
	private BranchEdge be[] = null; //the edges in the tree
	private int numEdges = 0; //the number of edges in the tree
	
	
	BranchTree(){
		bn = new BranchNode[50];
		numNodes = 0;
		
		be = new BranchEdge[50];
		numEdges = 0;
	}
	
	//Add a node to the tree
	public void addNode(BranchNode tn){
		
		if (numNodes>=bn.length) //the array must be extended
			bn = doubleSize(bn);
		
		bn[numNodes] = tn;
		bn[numNodes].setIndex(numNodes);
		
		numNodes++;
	}
	
	//Delete a node from the tree
	public void deleteNode(int index){
		
		deleteEdges(index); //delete all edges for the index node
		
		for (int i=index; i<(numNodes-1); i++){
			bn[i] = bn[i+1];
			bn[i].setIndex(i);
		}
		
		bn[numNodes-1] = null;
		numNodes--;
	}
	
	//Adds an edge between the nodes with indices i1 and i2
	public int addEdge(int i1, int i2){
		
		if (!edgeExists(i1,i2)) {//edge does not exist
			
			if (numEdges>=be.length) //array must be extended
				be = doubleSize(be);
			
			be[numEdges] = new BranchEdge(bn[i1],bn[i2]);
			bn[i1].addEdge();
			bn[i2].addEdge();
			numEdges++;
			
			return (numEdges-1);
		}
		return -1;
	}
	
	//Determines if an edge already exists between the nodes with indices i1 and i2
	public boolean edgeExists(int i1, int i2){
		
		for (int i=0; i<numEdges; i++){
			
			int n1 = be[i].getn1().getIndex();
			int n2 = be[i].getn2().getIndex();
			
			if ( ( (n1==i1) && (n2==i2) ) || ( (n1==i2) && (n2==i1) ) )
				return true;
		}
		return false;
	}
	
	//Delete an edge from the tree
	public void deleteEdge(int index){
		
		bn[be[index].getn1().getIndex()].deleteEdge();
		bn[be[index].getn2().getIndex()].deleteEdge();
		
		for (int i=index; i<(numEdges-1); i++){
			be[i] = be[i+1];
		}
		
		be[numEdges-1] = null;
		numEdges--;
	}
	
	//Delete all edges for the index'th node
	private void deleteEdges(int index){
		
		int i = 0;
		while (i<numEdges){
			
			if ( (be[i].getn1().getIndex()==index) || (be[i].getn2().getIndex()==index) ) //found an edge incident with bn
				deleteEdge(i);
			
			else
				i++;
		}
	}
	
	//Returns the index'th node from the node matrix
	public BranchNode getNode(int index){
		return bn[index];
	}
	
	public int getNumNodes(){
		return numNodes;
	}
	
	//Returns the index'th edge from the edge matrix
	public BranchEdge getEdge(int index){
		return be[index];
	}
	
	public int getNumEdges(){
		return numEdges;
	}
	
	//Returns the set of edges that involve the index'th node
	public BranchEdge [] getEdgesForNode(int index){
		
		BranchEdge bnEdges[] = new BranchEdge[bn[index].getNumEdges()];
		
		int c = 0;
		for (int i=0; i<numEdges; i++){
			
			if ( (be[i].getn1().getIndex()==index) || (be[i].getn2().getIndex()==index) ){
				bnEdges[c] = be[i];
				c++;
			}
			
			if (c>=bnEdges.length) //found all edges
				break;
		}
		
		return bnEdges;
	}
	
	//Returns the set of nodes adjacent to node bn in this graph;
	//If (directed==true), then only directed edges are counted; otherwise, undirected edges are counted
	public LinkedHashSet<BranchNode> getNeighborsForNode(BranchNode bn, boolean directed){
		
		LinkedHashSet<BranchNode> vN = new LinkedHashSet<BranchNode>();
		
		for (int i=0; i<numEdges; i++){
			
			if (be[i].getn1().getIndex()==bn.getIndex()) //n1 for current edge same as node bn
				vN.add(be[i].getn2()); //add the other node of the edge to the list of neighbors of bn
			
			else if ( (!directed) && (be[i].getn2().getIndex()==bn.getIndex()) ) //undirected, and n2 for current edges same as node bn
				vN.add(be[i].getn1()); //add the other node of the edge to the list of neighbors of bn
		}
		
		return vN;
	}
	
	//Double the size of the BranchEdge array
	private BranchEdge [] doubleSize(BranchEdge b[]){
		BranchEdge tmp[] = new BranchEdge[b.length*2];
		System.arraycopy(b, 0, tmp, 0, b.length);
		return tmp;
	}
	
	//Double the size of the BranchNode array
	private BranchNode [] doubleSize(BranchNode b[]){
		BranchNode tmp[] = new BranchNode[b.length*2];
		System.arraycopy(b, 0, tmp, 0, b.length);
		return tmp;
	}
	
	//Returns a deep copy of this BranchTree (modified from http://javatechniques.com/blog/faster-deep-copies-of-java-objects/, accessed 10/30/2008)
	public BranchTree deepCopy(){
		BranchTree c = null;
		try {
			ByteArrayOutputStream b = new ByteArrayOutputStream();
			ObjectOutputStream fout = new ObjectOutputStream(b);
			fout.writeObject(this);
			fout.flush();
			fout.close();
			
			ObjectInputStream fin = new ObjectInputStream(new ByteArrayInputStream(b.toByteArray()));
			c = (BranchTree)fin.readObject();
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}
		return c;
	}
}
