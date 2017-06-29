package edu.duke.cs.osprey.tools.branchdecomposition;
import java.io.*;

public class BranchNode implements Serializable {
	
	private boolean isLeaf = false; //leaf or internal node
	
	private String graphV1 = null; //first graph vertex associated with the current node
	private String graphV2 = null; //second graph vertex associated with the current node
	
	private int numEdges = 0; //number of edges incident with this node
	
	private int index = -1; //the index of this node into the node matrix
	
	
	BranchNode(boolean leafNode, String gv1, String gv2){
		
		isLeaf = leafNode;
		
		graphV1 = gv1;
		graphV2 = gv2;
	
		numEdges = 0;
		index = -1;
	}
	
	public String getv1(){
		return graphV1;
	}
	
	public String getv2(){
		return graphV2;
	}
	
	public int getNumEdges(){
		return numEdges;
	}
	
	public boolean getIsLeaf(){
		return isLeaf;
	}
	
	public int getIndex(){
		return index;
	}
	
	public void setIndex(int i){
		index = i;
	}
	
	//Updates the number of edges (adds an edge) for this node
	public void addEdge(){		
		numEdges++;
	}
	
	//Updates the number of edges (deletes and edge) for this node
	public void deleteEdge(){
		numEdges--;
	}
}
