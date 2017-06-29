package edu.duke.cs.osprey.tools.branchdecomposition;
import java.util.*;
import java.io.*;

public class BranchEdge implements Serializable {

	private BranchNode n1 = null; //the first node adjacent to the edge in the rooted tree
	private BranchNode n2 = null; //the second node adjacent to the edge in the rooted tree
	
	private LinkedHashSet<String> M = null; //the M set
	
	private boolean checkedn1 = false; //determines if this edge has been checked for pushing toward node n1
	private boolean checkedn2 = false; //determines if this edge has been checked for pushing toward node n2
	
	BranchEdge(BranchNode bn1, BranchNode bn2){
		
		n1 = bn1;
		n2 = bn2;
		
		M = new LinkedHashSet<String>();
		
		checkedn1 = false;
		checkedn2 = false;
	}
	
	public LinkedHashSet<String> getM(){
		return M;
	}
	
	public BranchNode getn1(){
		return n1;
	}
	
	public BranchNode getn2(){
		return n2;
	}
	
	public void setM(LinkedHashSet<String> nM){
		M = new LinkedHashSet<String>(nM);
	}
	
	//Add v to M
	public void addM(String v){
		M.add(v);
	}
	
	//Determines if this edge has been checked for pushing toward node bn 
	public boolean isChecked(BranchNode bn){
		
		if (bn.getIndex()==n1.getIndex())
			return checkedn1;
		
		else if (bn.getIndex()==n2.getIndex())
			return checkedn2;
		
		else {
			System.out.println("ERROR: current edge not incident to node "+bn.getIndex());
			System.exit(1);
			return false;
		}
	}
	
	//Sets the checked flag for node bn
	public void setChecked(BranchNode bn, boolean c){
		
		if (bn.getIndex()==n1.getIndex())
			checkedn1 = c;
		
		else if (bn.getIndex()==n2.getIndex())
			checkedn2 = c;
		
		else {
			System.out.println("ERROR: current edge not incident to node "+bn.getIndex());
			System.exit(1);
		}
	}
}
