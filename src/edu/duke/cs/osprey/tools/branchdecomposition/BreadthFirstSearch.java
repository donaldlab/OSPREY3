package edu.duke.cs.osprey.tools.branchdecomposition;
import java.util.*;
import java.util.concurrent.*;

public class BreadthFirstSearch {
	
	private final String white = "w";
	private final String gray = "g";
	private final String black = "b";

	private String color[] = null; //the color of each node (w, g, or b)	
	private int d[] = null; //the depth of each node
	private int pi[] = null; //the parent (in the search) of each node
	
	private BranchTree G = null; //the directed graph to be searched
	private int s = -1; //the index of the source node
	private int t = -1; //the index of the end node
	private int adj[][] = null; //the list of adjacent nodes for each node (the direction of edges is taken into account)
	
	private ArrayBlockingQueue<Integer> q = null; //the current queue
	
	BreadthFirstSearch(BranchTree graph, int sInd, int tInd){
		
		G = graph;
		s = sInd;
		t = tInd;
		int numV = G.getNumNodes();
		
		adj = new int[numV][];
		int curInd[] = new int[numV];
		for (int i=0; i<G.getNumEdges(); i++){ //compute the adjacency lists and the vertices-edges map
			
			BranchEdge curE = G.getEdge(i);
			int n1 = curE.getn1().getIndex();
			int n2 = curE.getn2().getIndex();
			
			if (adj[n1]==null){
				adj[n1] = new int[curE.getn1().getNumEdges()];
				curInd[n1] = 0;
			}
			adj[n1][curInd[n1]] = n2;
			curInd[n1]++;
			
			if (adj[n2]==null){
				adj[n2] = new int[curE.getn2().getNumEdges()];
				curInd[n2] = 0;
			}
			adj[n2][curInd[n2]] = n1;
			curInd[n2]++;
		}
	}
	
	//Returns the next augmenting path from node s to node t by also checking the residual capacities of the edges
	public int [] getNextPath(int rc[][]){
		
		init();
		
		while (!q.isEmpty()){ //more elements in the queue
			
			int u = q.remove();
			for (int vc=0; vc<adj[u].length; vc++){
				
				int v = adj[u][vc];
				
				if (v==t) { //the end node, so return the path
					pi[v] = u;
					return tracePath(u);
				}
				
				else if (color[v].equalsIgnoreCase(white)){ //node v not discovered yet
					
					if (rc[u][v]>0) { //there is residual capacity for this edge
						
						color[v] = gray;
						d[v] = d[u] + 1;
						pi[v] = u;
						q.add(v);
					}
				}
			}
			color[u] = black;
		}
		return null;
	}
	
	//Traces the path involvind the end node t and its parent (in the current state of the search) pOft, back to the source node s;
	//		Must be called only when a full path is computed in getNextPath()
	private int [] tracePath(int pOft){
		
		LinkedHashSet<Integer> p = new LinkedHashSet<Integer>();
		
		int curNode = t; //start with the end node
		while (curNode!=s){ //not reached the source node
			p.add(curNode); //add to path
			curNode = pi[curNode]; //move to the parent of curNode
		}
		p.add(s);
		
		Object pA[] = p.toArray(); //this is the reverse path
		
		int pAr[] = new int [pA.length];
		for (int i=0; i<pAr.length; i++){ //find the forward path
			pAr[i] = (Integer)pA[pA.length-1-i];
		}
		
		return pAr;
	}
	
	//Initializes color[], d[], and pi[], and the queue q
	public void init(){
		
		int numV = G.getNumNodes();
		
		color = new String[numV];
		d = new int[numV];
		pi = new int[numV];
		
		for (int i=0; i<color.length; i++){ //initialize
			color[i] = white;
			d[i] = (int)Math.pow(10, 10); //this should be greater than any depth in the search
			pi[i] = -1;
		}
		color[s] = gray;
		
		q = new ArrayBlockingQueue<Integer>(numV);
		q.add(s);
	}
}
