package edu.duke.cs.osprey.tools.branchdecomposition;
import java.util.*;

/*
 * The vertex information for the original graph
 */
public class GraphVertices {
	
	private String vNames[] = null; //the names of the vertices
	int numV = 0; //the number of vertices
	
	private LinkedHashSet<String> graphVertices = null; //the set of graph vertices

	GraphVertices(){		
		vNames = new String[10];
		numV = 0;
		graphVertices = null;
	}
	
	//Adds the vertex name to the graph;
	//Returns false if this vertex had already been added to the graph, and true otherwise
	public boolean addV(String name){
		
		if (!isAlreadyAdded(name)){ //vertex not already added, so add
			
			if (numV>=vNames.length){ //increase the size of the array
				
				String tmps[] = new String[vNames.length*2];
				System.arraycopy(vNames, 0, tmps, 0, vNames.length);
				vNames = tmps;
			}
			
			vNames[numV] = name;
			numV++;
			
			return true;
		}
		
		else
			return false;
	}
	
	//Checks if the vertex specified by name has already been added to the graph
	private boolean isAlreadyAdded(String name){
		for (int i=0; i<numV; i++){
			if (vNames[i].equalsIgnoreCase(name)) //vertex already in graph
				return true;
		}
		return false;
	}
	
	//The graph is complete (no more vertices will be added), so set the LinkedHashSet graphVertices
	public void setCompleteGraph(){
		graphVertices = new LinkedHashSet<String>();
		for (int i=0; i<numV; i++)
			graphVertices.add(vNames[i]);
	}
	
	public LinkedHashSet<String> getGraphVertices(){
		return graphVertices;
	}
	
	public int getNumV(){
		return numV;
	}
	
	//Returns the index into vNames[] for vertex name
	public int getVind(String name){
		
		for (int i=0; i<numV; i++){
			if (vNames[i].equalsIgnoreCase(name)) //found vertex
				return i;
		}
		
		System.out.println("ERROR: graph vertex "+name+" not found");
		System.exit(1);
		return -1;
	}
}
