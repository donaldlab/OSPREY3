package edu.duke.cs.osprey.sparse;

import java.io.File;

public class TreeNode {
	private int name = -1; //the name of this node
	
	private boolean isLeaf = false; //leaf or internal node
	
	private int graphV1 = -1; //first graph vertex associated with the current node (molecule-relative numbering)
	private int graphV2 = -1; //second graph vertex associated with the current node (molecule-relative numbering)

	private TreeNode p = null; //parent node in the tree
	private TreeNode lc = null; //left child of this node
	private TreeNode rc = null; //right child of this node
	
	private TreeEdge cOfEdge = null; //child of which tree edge is this node (the other end of this edge is adjacent to the parent 'p' of this node)


	public TreeNode(int tName, boolean tIsLeaf, int v1, int v2){
		
		name = tName;
		isLeaf = tIsLeaf;
		graphV1 = v1;
		graphV2 = v2;
		
		p = null;
		lc = null;
		rc = null;
		cOfEdge = null;
	}
	
	public TreeNode getp(){
		return p;
	}
	
	public TreeNode getlc(){
		return lc;
	}
	
	public TreeNode getrc(){
		return rc;
	}
	
	public int getv1(){
		return graphV1;
	}
	
	public int getv2(){
		return graphV2;
	}
	
	public TreeEdge getCofEdge(){
		return cOfEdge;
	}
	
	public void setCofEdge(TreeEdge ce){
		cOfEdge = ce;
	}
	
	public void setLc(TreeNode slc){
		lc = slc;
	}
	
	public void setRc(TreeNode src){
		rc = src;
	}
	
	public void setP(TreeNode sp){
		p = sp;
	}
	
	public int getName(){
		return name;
	}
	
	public boolean isLeaf(){
		return isLeaf;
	}

	public String toString()
	{
		String out = "[";
		if(cOfEdge != null)
		{			
			for(int i : cOfEdge.getLambda())
				out+=i+", ";
		}
		out+="]";
		/*
		out += " - L Set:[";
		if(cOfEdge != null)
		{			
			for(int i : cOfEdge.getL())
				out+=i+", ";
		}
		out+="]";
		
		out += " - M Set:[";
		if(cOfEdge != null)
		{			
			for(int i : cOfEdge.getM())
				out+=i+", ";
		}
		
		out+="]";
		*/
		return out + "id:" + this.hashCode();
	}
	
    public void printTree(String prefix)
    {
        String output = prefix+this;

        System.out.println(output);
        
        if(cOfEdge != null && cOfEdge.getLambda().containsAll(cOfEdge.getL()))
        	return;

        
        if(lc != null)
        lc.printTree(prefix+"L--");
        if(rc!=null)
        rc.printTree(prefix+"R--");
    }

}
