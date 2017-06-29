package edu.duke.cs.osprey.tools.branchdecomposition;
import java.util.*;

public class MinVertexCut {
	
	private LinkedHashSet<String> Lset = null;
	private LinkedHashSet<String> cutSet = null;
	private LinkedHashSet<String> Rset = null;
	
	MinVertexCut(){
	}
	
	//Find the min vertex cut for undirected graph H, from source node si to sink node ti (indexed in the node array of H);
	//		The cut set is stored in the set cutSet, and the resulting two partitions of the remaining vertices are stored in Lset and Rset;
	//The graph H is not modified here
	public void findMinVertexCut(BranchTree H, int si, int ti){
		
		Lset = new LinkedHashSet<String>(); 
		cutSet = new LinkedHashSet<String>(); 
		Rset = new LinkedHashSet<String>();
		
		//First, transform H into a directed graph suitable for the max flow algorithm
		int vMap[][] = new int[H.getNumNodes()][2];
		int splitV[] = new int[2*H.getNumNodes()];
		BranchTree dH = transformGraph(H.deepCopy(),vMap,splitV);
		
		//The indices of the source and the sink in the directed graph dH
		int siNew = vMap[si][1];
		int tiNew = vMap[ti][0];
		
		//Compute the min cut (using max flow) of dH and return the corresponding vertex partitions into the sets S1 and S2
		LinkedHashSet<Integer> S1 = new LinkedHashSet<Integer>();
		LinkedHashSet<Integer> S2 = new LinkedHashSet<Integer>();
		compMinCut(dH, siNew, tiNew, S1, S2, splitV);
		
		//Compute the minimum vertex cut and the corresponding partition of the vertices into three sets;
		//		The minimum vertex cut of H consists of all vertices v in the original graph H, such that the corresponding
		//		vertices v' and v" from dH are not both in S1 and not both in S2;
		//Always add the source to Lset and the sink to Rset
		for (int i=0; i<H.getNumNodes(); i++){	
			
			if (i==si) //the source
				Lset.add(H.getNode(i).getv1());
			
			else if (i==ti) //the sink
				Rset.add(H.getNode(i).getv1());
			
			else { //not the source or the sink
				if ( S1.contains(vMap[i][0]) && S1.contains(vMap[i][1])) //both dH vertices are in S1
					Lset.add(H.getNode(i).getv1());
				
				else if ( S2.contains(vMap[i][0]) && S2.contains(vMap[i][1])) //both dH vertices are in S2
					Rset.add(H.getNode(i).getv1());
				
				else //one vertex is in S1 and the other one is in S2
					cutSet.add(H.getNode(i).getv1());
			}
		}
		
		Lset.remove(""); //remove entries not corresponding to graph vertices
		cutSet.remove("");
		Rset.remove("");
	}
	
	//Transform the undirected graph H into a directed graph dH;
	//		For each vertex v in H, two vertices v' and v" and a directed edge (v',v") are added in dH;
	//		For each undirected edge (u,v) in H, two directed edges (u",v') and (v",u') are added in dH;
	//		The mapping between the indices of each vertex of H and the two corresponding vertices in the new graph is returned in vMap[][];
	//		The mapping between v' and v" (and, respectively, v" and v') is returned in splitV[]
	private BranchTree transformGraph(BranchTree H, int vMap[][], int splitV[]){
		
		BranchTree dH = new BranchTree();
		
		for (int i=0; i<H.getNumNodes(); i++){ //add nodes based on H
			
			BranchNode nodeH = H.getNode(i);
			
			BranchNode bn1 = new BranchNode(nodeH.getIsLeaf(),nodeH.getv1(),nodeH.getv2());
			dH.addNode(bn1); //add the first node
			vMap[i][0] = bn1.getIndex();
			
			BranchNode bn2 = new BranchNode(nodeH.getIsLeaf(),nodeH.getv1(),nodeH.getv2());
			dH.addNode(bn2); //add the second node
			vMap[i][1] = bn2.getIndex();
			
			splitV[vMap[i][0]] = vMap[i][1];
			splitV[vMap[i][1]] = vMap[i][0];
			
			dH.addEdge(vMap[i][0], vMap[i][1]); //add the edge from node bn1 to bn2 (this is a directed edge)
		}
		
		for (int j=0; j<H.getNumEdges(); j++){ //add the new directed edges based on the undirected edges in H
			
			BranchEdge be = H.getEdge(j);
			int i1 = be.getn1().getIndex();
			int i2 = be.getn2().getIndex();
			
			dH.addEdge(vMap[i1][1], vMap[i2][0]);
			dH.addEdge(vMap[i2][1], vMap[i1][0]);
		}
		
		return dH;
	}
	
	//Computes the min cut from the source node with index s to the sink node with index t in the directed graph dH;
	//		Outputs the cut (a partition of the vertices) in the sets S1 and S2
	private void compMinCut(BranchTree dH, int s, int t, LinkedHashSet<Integer> S1, LinkedHashSet<Integer> S2, int splitV[]){
		
		//Assign unit capacity to the edges in dH, and initialize the flow to 0;
		int c[][] = new int[dH.getNumNodes()][dH.getNumNodes()]; //initial capacities
		int f[][] = new int[c.length][c.length];
		int rc[][] = new int[c.length][c.length]; //residual capacities
		
		for (int i=0; i<c.length; i++){
			for (int j=0; j<c[i].length; j++){
				c[i][j] = 0;
				f[i][j] = 0;
				rc[i][j] = c[i][j];
			}
		}
		
		for (int i=0; i<dH.getNumEdges(); i++){ //set the capacities for all (directed) edges in dH to 1
			
			BranchEdge curE = dH.getEdge(i);
			int n1 = curE.getn1().getIndex();
			int n2 = curE.getn2().getIndex();
			
			c[n1][n2] = 1;
			rc[n1][n2] = c[n1][n2]; //set the residual capacity to the initial capacity
		}
		
		int p[] = null; //used to store the path output from the search
		BreadthFirstSearch bfs = new BreadthFirstSearch(dH,s,t);
		while ( (p=bfs.getNextPath(rc)) != null ){ //there are more paths
			
			int rcp = 10000000; //the residual capacity of path p
			for (int i=0; i<p.length-1; i++)
				rcp = Math.min(rcp, rc[p[i]][p[i+1]]);
			
			for (int i=0; i<p.length-1; i++){
				
				f[p[i]][p[i+1]] += rcp; //flow for this edge
				f[p[i+1]][p[i]] = -f[p[i]][p[i+1]]; //flow for the opposite edge
			}
			
			for (int i=0; i<p.length-1; i++){ //update the residual capacities
				rc[p[i]][p[i+1]] = c[p[i]][p[i+1]] - f[p[i]][p[i+1]];
				rc[p[i+1]][p[i]] = c[p[i+1]][p[i]] - f[p[i+1]][p[i]];
			}
		}
		
		
		//Find all nodes in the layered network for the resulting (max) flow
		int lnNodes[] = getLayeredNetworkNodes(dH,s,t,c,f,splitV);
		for (int i=0; i<lnNodes.length; i++)
			S1.add(lnNodes[i]);
		
		for (int i=0; i<dH.getNumNodes(); i++)
			S2.add(i);
		S2.removeAll(S1); //the second partition
	}
	
	//Constructs the layered network corresponding to the graph dH and the current flow in f[][], and returns all nodes in that layered network
	private int [] getLayeredNetworkNodes(BranchTree dH, int s, int t, int c[][], int f[][], int splitV[]){
		
		for (int i=0; i<dH.getNumEdges(); i++){ //set the capacities for all external edges in dH to infinity
			
			BranchEdge curE = dH.getEdge(i);
			int n1 = curE.getn1().getIndex();
			int n2 = curE.getn2().getIndex();
			
			if (splitV[n1]!=n2)  //external edge in dH, set its capacity to infinity
				c[n1][n2] = Integer.MAX_VALUE;		
		}

		//Construct the layers
		int Vi[][] = new int[dH.getNumNodes()][];
		Vi[0] = new int[1];
		Vi[0][0] = s;
		
		int Vall[] = new int[dH.getNumNodes()];
		int curAll = 0;
		Vall[curAll] = Vi[0][0];
		curAll++;
		
		int ind = 0;
		boolean done = false;
		while (!done){ //sink t not reached yet
			
			int T[] = new int[dH.getNumNodes()];
			int curT = 0;
			
			for (int i=0; i<Vi[ind].length; i++){
				
				for (int j=0; j<dH.getNumNodes(); j++){
					
					if (!isElement(j,Vall,curAll)){ //j is not part of any layer yet
						
						if ( (f[Vi[ind][i]][j]<c[Vi[ind][i]][j]) || (f[j][Vi[ind][i]]>0) ) { //the edge between the current node of Vi and node j of dH is useful
							
							T[curT] = j;
							curT++;
							
							Vall[curAll] = j;
							curAll++;
						}
					}
				}
			}
			
			if (curT==0) //no more layers, so done
				done = true;
			
			else if (isElement(t,T,curT)){ //the sink t is an element of this layer, update the next layer, and done
				
				Vi[ind+1] = new int[1];
				Vi[ind+1][0] = t;
				
				done = true;
			}
			
			else { //update and move to the next layer
				int tmp[] = new int[curT];
				System.arraycopy(T, 0, tmp, 0, curT);
				T = tmp;
				
				Vi[ind+1] = T;
				ind++;
			}
		}
		
		int tmp[] = new int[curAll];
		System.arraycopy(Vall, 0, tmp, 0, curAll);
		Vall = tmp;
		
		return Vall; //return all nodes in the layered network
	}
	
	//Checks if a is found among the first sizeA elements of A[]
	private boolean isElement(int a, int A[], int sizeA){
		for (int i=0; i<sizeA; i++){
			if (A[i]==a)
				return true;
		}
		return false;
	}
	
	public LinkedHashSet<String> getCutSet(){
		return cutSet;
	}
	
	public LinkedHashSet<String> getLset(){
		return Lset;
	}
	
	public LinkedHashSet<String> getRset(){
		return Rset;
	}
	
	public void addToLset(String s){
		Lset.add(s);
	}
	
	public void addToRset(String s){
		Rset.add(s);
	}
}
