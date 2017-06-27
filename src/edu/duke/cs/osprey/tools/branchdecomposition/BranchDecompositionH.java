package edu.duke.cs.osprey.tools.branchdecomposition;
import java.io.*;
import java.util.*;
import java.lang.management.*;
//import static org.math.array.LinearAlgebra.*;
import Jama.*;


public class BranchDecompositionH {
	
	private static final boolean debug = false;
	
	private BranchTree bt = null; //the branch decomposition tree
	
	GraphVertices gv = null; //the set of graph vertices
	
	public BranchDecompositionH(String args[]){
		
		long startTimeCPU = CPUTime();
		long startTimeUser = UserTime();
		long startTimeSystem = SystemTime();
		
		LinkedHashSet<String> gvDup = new LinkedHashSet<String>(); //the set of graph vertices found in more than one node
		readGraphFile(args[0],gvDup); //read in the residue interaction graph
		
		constructStarGraph(gvDup); //construct the initial star graph
		
		boolean done = false;
		boolean done2sep = false;
		boolean done3sep = false;
		while (!done){ //perform a new split each iteration until a branch decomposition is obtained
			
			LinkedHashSet<BranchNode> bigNodes = getNodesToSplit(); //the set of nodes to split (all nodes with degree greater than three) 
			
			if (!bigNodes.isEmpty()){ //only partial branch decomposition (more nodes with degree greater than three)
				
				System.out.println("Starting a new split.. nodes with degree >3: "+bigNodes.size());
				
				LinkedHashSet<BranchEdge> sx = new LinkedHashSet<BranchEdge>();
				LinkedHashSet<BranchEdge> sy = new LinkedHashSet<BranchEdge>();
				
				BranchNode pn = findPushNode(bigNodes); //get a node to push (if any such nodes exist)
				
				if (pn!=null){ //try pushing on node pn, since it has edges that have not been checked for pushing yet
					
					System.out.println("Attempting to push on node "+pn.getIndex()+" ("+pn.getNumEdges()+" edges)");
					
					pushNode(pn,sx,sy);
				}
				
				else if (!done2sep){ //no pushing candidate found, check for 2-separations
					
					boolean found = false;
					
					Object bigNodesA[] = bigNodes.toArray();
					for (int i=0; i<bigNodesA.length; i++){ //check all nodes with degree>3 for 2-separations
						
						System.out.println("Pushing not possible.. attempting 2-separation at node "+((BranchNode)bigNodesA[i]).getIndex()+" ("+((BranchNode)bigNodesA[i]).getNumEdges()+" edges)");
						
						if (find2sep((BranchNode)bigNodesA[i],sx,sy)){ //2-separation found
							pn = (BranchNode)bigNodesA[i];
							found = true;
							break;
						}
					}
					
					if (!found) //no more 2-separations exist, so no further checking for 2-separations is necessary
						done2sep = true;
				}
				
				else if (!done3sep){ //no pushing candidate found and no more 2-separations possible, so check for 3-separations
					
					boolean found = false;
					
					Object bigNodesA[] = bigNodes.toArray();
					for (int i=0; i<bigNodesA.length; i++){ //check all nodes with degree>3 for 3-separations
						
						System.out.println("Pushing not possible.. attempting 3-separation at node "+((BranchNode)bigNodesA[i]).getIndex()+" ("+((BranchNode)bigNodesA[i]).getNumEdges()+" edges)");
						
						if (find3sep((BranchNode)bigNodesA[i],sx,sy)) { //3-separation found
							pn = (BranchNode)bigNodesA[i];
							found = true;
							break;
						}
					}
					
					if (!found) //no more 3-separations exist, so no further checking for 3-separations is necessary
						done3sep = true;
				}
				
				else { //no more 2- or 3-separations, so perform the eigenvector heuristic
					
					pn = (BranchNode)bigNodes.toArray()[0]; //find a node with degree greater than three for splitting 
					
					System.out.println("Pushing not possible.. performing eigenvalue heuristic on node "+pn.getIndex()+" ("+pn.getNumEdges()+" edges)");
					
					doEigen(pn,sx,sy);
				}
				
				if (!sx.isEmpty()) { //successfully found a partition, so split
					
					System.out.print("Partitioning successful.. performing split..");
					
					performSplit(pn,sx,sy);
					
					if (debug)
						outputBranchDecomposition(args[1]+"_"+System.currentTimeMillis());
				}
			}
			
			else { //branch decomposition obtained, so done
				done = true;
			}
		}
		
		outputBranchDecomposition(args[1]);
		
		long endTimeCPU = CPUTime();
		long endTimeUser = UserTime();
		long endTimeSystem = SystemTime();
		System.out.println("DONE..");
		System.out.println("Time: CPU "+((endTimeCPU-startTimeCPU)/60000000000.0));
		System.out.println("Time: User "+((endTimeUser-startTimeUser)/60000000000.0));
		System.out.println("Time: System "+((endTimeSystem-startTimeSystem)/60000000000.0));
	}
	
	//Construct the initial star graph; must be called immediately after reading in the residue interaction graph
	private void constructStarGraph(LinkedHashSet<String> gvDup){
		
		BranchNode starNode = new BranchNode(false,null,null);
		bt.addNode(starNode); //create the central node
		int starNodeIndex = starNode.getIndex();
		
		for (int i=0; i<bt.getNumNodes(); i++){ //create an edge between the central node and every other node
			
			if (i!=starNodeIndex){
				
				int eInd = bt.addEdge(i,starNodeIndex); //add the edge
				
				
				//Determine the M, L, and R sets for the new edge
				BranchEdge be = bt.getEdge(eInd);
				
				BranchNode bn = bt.getNode(i);
				String gv1 = bn.getv1();
				String gv2 = bn.getv2();
				
				LinkedHashSet<String> curR = new LinkedHashSet<String>(gv.getGraphVertices()); //the R set contains every vertex other than gv1 and gv2
				curR.remove(gv1);
				curR.remove(gv2);
				
				LinkedHashSet<String> curL = new LinkedHashSet<String>(); //the L set may contain none, one, or both of gv1 and gv2
				curL.add(gv1);
				curL.add(gv2);
				
				if (gvDup.contains(gv1)){ //gv1 should be in M
					be.addM(gv1);
					curL.remove(gv1);
				}
				if (gvDup.contains(gv2)){ //gv2 should be in M
					be.addM(gv2);
					curL.remove(gv2);
				}
			}
		}
	}
	
	//Returns all nodes (and the adjacent edges that have not been checked for pushing) with degree greater than three
	private LinkedHashSet<BranchNode> getNodesToSplit(){
		
		LinkedHashSet<BranchNode> pn = new LinkedHashSet<BranchNode>();
		
		for (int i=0; i<bt.getNumNodes(); i++){
			BranchNode bn = bt.getNode(i);			
			if ( (!bn.getIsLeaf()) && (bn.getNumEdges()>3) ) //internal node with degree greater than three				
				pn.add(bn);				
		}
		return pn;
	}
	
	//Determine if any edges for any of the nodes in the set bigNodes have not already been checked for pushing;
	//		Return one such node, or null if no such node exists
	private BranchNode findPushNode(LinkedHashSet<BranchNode> bigNodes){
		
		Iterator<BranchNode> bigNodesIt = bigNodes.iterator();
		while (bigNodesIt.hasNext()){
			BranchNode cn = bigNodesIt.next();
			BranchEdge ce[] = bt.getEdgesForNode(cn.getIndex());
			for (int j=0; j<cn.getNumEdges(); j++){
				if (!ce[j].isChecked(cn))
					return cn; //found one such node, so return it
			}
		}
		return null; //no such nodes (edges) exist
	}
	
	//Try pushing on node pn;
	//		if a pair of edges is found to satisfy the pushing inequality, perform and return the edge partition in the sets sx and sy
	private void pushNode(BranchNode pn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy){
		
		BranchEdge pe[] = bt.getEdgesForNode(pn.getIndex()); //all the edges adjacent to node pn
		
		boolean vInMe[][] = new boolean[gv.getNumV()][pe.length]; //the edges from pe[] for which each vertex belongs to the M set
		for (int k=0; k<pe.length; k++){
			
			Object curMe[] = pe[k].getM().toArray();
			for (int i=0; i<curMe.length; i++){ //update vInMe[][] for all the vertices in curMe[]
				vInMe[gv.getVind((String)curMe[i])][k] = true;
			}
		}
		
		for (int i=0; i<pe.length; i++){
			if (!pe[i].isChecked(pn)){ //an unchecked edge
				
				pe[i].setChecked(pn,true); //set the checked flag
				
				for (int j=0; j<pe.length; j++){
					if (i!=j){
						
						LinkedHashSet<String> uMij = new LinkedHashSet<String>(pe[i].getM());
						uMij.addAll(pe[j].getM()); //the union of the M sets for edges i and j	
						Object uMijA[] = uMij.toArray();
						
						//Delete from uMij all vertices that are found only in the union of the M sets for edges i and j,
						//		but are not found in the M sets for any other edges in pe[]
						for (int k=0; k<uMijA.length; k++){
							
							int curVind = gv.getVind((String)uMijA[k]);							
							boolean found = false;
							for (int l=0; l<vInMe[curVind].length; l++){
								if ( (l!=i) && (l!=j) && (vInMe[curVind][l]) ){ //the current vertex is in the M set for an edge other than edges i and j
									found = true;
									break;
								}
							}
							
							if (!found) { //the current vertex is only in the M sets of edges i and/or j, but in no other edges in pe[], so delete
								uMij.remove((String)uMijA[k]);
							}
						}
						
						//At this point, uMij represents the set of vertices in the intersection of two unions:
						//		the union of the M sets of edges i and j, and
						//		the union of the M sets of all edges in pe[] other than i and j;
						//Check the pushing inequality
						if ( uMij.size() <= Math.max(pe[i].getM().size(), pe[j].getM().size()) ){ //pushing inequality holds, so perform the partition and return
							
							sx.add(pe[i]);
							sx.add(pe[j]);
							
							for (int k=0; k<pe.length; k++){
								if ( (k!=i) && (k!=j) )
									sy.add(pe[k]);
							}
							
							return;
						}
					}
				}
			}
		}
	}
	
	//Tries to find a 2-separation for node pn; the resulting partition is returned in the sets sx and sy;
	//	Returns true if a 2-separation is found, and false otherwise
	private boolean find2sep(BranchNode pn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy){
		
		BranchEdge D[] = bt.getEdgesForNode(pn.getIndex()); //all the edges incident with node pn
		
		LinkedHashSet<String> S = new LinkedHashSet<String>(); //the union of the M sets for all edges in D[]
		for (int i=0; i<D.length; i++){
			S.addAll(D[i].getM());
		}		
		Object Sa[] = S.toArray();
		
		//Construct a new graph H with nodes the vertices in the union of the M sets for all edges in D[],
		//		and an edge between a pair of nodes if some M set for an edge in D[] contains both these nodes
		BranchTree H = constructGraphH(Sa,D);
		
		return find23sepHelper(H,sx,sy,D,true);
	}
	
	//Tries to find a 3-separation for node pn; the resulting partition is returned in the sets sx and sy;
	//	Returns true if a 3-separation is found, and false otherwise
	private boolean find3sep(BranchNode pn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy){
		
		BranchEdge D[] = bt.getEdgesForNode(pn.getIndex()); //all the edges incident with node pn
		
		LinkedHashSet<String> S = new LinkedHashSet<String>(); //the union of the M sets for all edges in D[]
		for (int i=0; i<D.length; i++){
			S.addAll(D[i].getM());
		}		
		Object Sa[] = S.toArray();
		
		//Construct a new graph H with nodes the vertices in the union of the M sets for all edges in D[],
		//		and an edge between a pair of nodes if some M set for an edge in D[] contains both these nodes
		BranchTree H = constructGraphH(Sa,D);
		
		
		//Find a node x1 in H, such that x1 is of degree three and is part of a triangle; if such a node exists, perform the split
		for (int i=0; i<H.getNumNodes(); i++){
			BranchNode x1 = H.getNode(i);
			if (x1.getNumEdges()==3) { //x1 is of degree three
				for (int j=0; j<H.getNumNodes(); j++){ //for all nodes adjacent to x1
					if (j!=i) {
						BranchNode x2 = H.getNode(j);
						if (H.edgeExists(x1.getIndex(),x2.getIndex())){ //there is an edge between x1 and x2
							for (int k=j+1; k<H.getNumNodes(); k++){ //look for a triangle
								if (k!=i) {
									BranchNode x3 = H.getNode(k);
									if ( (H.edgeExists(x1.getIndex(),x3.getIndex())) && (H.edgeExists(x2.getIndex(),x3.getIndex())) ) { //found a triangle
										LinkedHashSet<BranchEdge> sxtmp = new LinkedHashSet<BranchEdge>();
										LinkedHashSet<BranchEdge> sytmp = new LinkedHashSet<BranchEdge>();
										for (int l=0; l<D.length; l++){
											sytmp.add(D[l]);
											if (D[l].getM().contains(x1.getv1())) //x1 is in the M set of the current edge, so add to the set sxtmp
												sxtmp.add(D[l]);
										}
										sytmp.removeAll(sxtmp); //sxtmp and sytmp are now a partition of D[]
										if (sytmp.size()>=2){ //a valid partition, so assign sx and sy and return
											sx.addAll(sxtmp);
											sy.addAll(sytmp);
											System.out.println("3-sp triangle: sx("+sx.size()+"), sy("+sy.size()+")");
											return true;
										}
										else { //not a valid partition
											System.out.println("3-sp triangle: cannot partition, size of sy is "+sy.size());
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		
		//Check for 3-separations with the L and R sets of size at least two
		return find23sepHelper(H,sx,sy,D,false);
	}
	
	//Computes the min vertex connectivity for graph H and performs a check if there is a 3-separation (such that the L and R sets each have at least two elements),
	//		or if there is a 2-separation (do2sep==true)
	private boolean find23sepHelper(BranchTree H, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy, BranchEdge D[], boolean do2sep) {
		
		//Find a node in H with minimum degree
		BranchNode nn = null;
		int minDegree = Integer.MAX_VALUE;
		for (int i=0; i<H.getNumNodes(); i++){
			if (H.getNode(i).getNumEdges()<minDegree){
				nn = H.getNode(i);
				minDegree = H.getNode(i).getNumEdges();
			}
		}
		
		System.out.println("min vertex degree: "+minDegree);
		
		//Find all neighbors of node nn (all nodes adjacent to node nn in graph H)
		BranchEdge cDedges[] = H.getEdgesForNode(nn.getIndex()); //the edges in H incident with node nn
		BranchNode nnNeighborsN[] = new BranchNode[cDedges.length];
		LinkedHashSet<Integer> nnNeighbors = new LinkedHashSet<Integer>(); //the vertices in H adjacent to node nn
		for (int i=0; i<cDedges.length; i++){
			if (cDedges[i].getn1().getIndex()==nn.getIndex()){
				nnNeighbors.add(cDedges[i].getn2().getIndex());
				nnNeighborsN[i] = cDedges[i].getn2();
			}
			else {
				nnNeighbors.add(cDedges[i].getn1().getIndex());
				nnNeighborsN[i] = cDedges[i].getn1();
			}
		}
		
		//For all non-neighbors of node nn in graph H, compute the vertex connectivity and check if it is a 2-separation
		int vConn = Integer.MAX_VALUE;
		for (int i=0; i<H.getNumNodes(); i++){
			if ( (!nnNeighbors.contains(H.getNode(i).getIndex())) && (H.getNode(i).getIndex()!=nn.getIndex()) ) { //for each non-neighbor of node nn
				
				BranchEdge e = new BranchEdge(new BranchNode(false,null,null), new BranchNode(false,null,null));
				
				if ( (do2sep) && (find2sepHelper(H,nn,H.getNode(i),sx,sy,D,e)) ) //found a 2-separation
					return true;
				
				if ( (!do2sep) && (findKsepPairV(H,nn,H.getNode(i),3,sx,sy,D,e)) ) //found a 3-separation
					return true;
				
				if (e.getM().size()>0)
					vConn = Math.min(vConn, e.getM().size());
			}
		}
		
		for (int i=0; i<nnNeighborsN.length; i++) {
			for (int j=i+1; j<nnNeighborsN.length; j++) {
				if (!H.edgeExists(nnNeighborsN[i].getIndex(), nnNeighborsN[j].getIndex())) { //no edge between nodes i and j in H, so compute the local vertex connectivity
					
					BranchEdge e = new BranchEdge(new BranchNode(false,null,null), new BranchNode(false,null,null));
					
					if ( (do2sep) && (find2sepHelper(H,nnNeighborsN[i],nnNeighborsN[j],sx,sy,D,e)) ) //found a 2-separation
						return true;
					
					if ( (!do2sep) && (findKsepPairV(H,nnNeighborsN[i],nnNeighborsN[j],3,sx,sy,D,e)) ) //found a 3-separation
						return true;
					
					if (e.getM().size()>0)
						vConn = Math.min(vConn, e.getM().size()); //update the min vertex connectivity
				}
			}
		}
		
		System.out.println("min vertex connectivity: "+vConn);
		
		return false;
	}
	
	//Helper for find2sep(); computes and returns (in edge e) the M set for a separation of graph H;
	//		The sets sxtmpV and sytmpV, and the graph H are not modified here
	private boolean find2sepHelper(BranchTree H, BranchNode v, BranchNode w, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy, BranchEdge D[], BranchEdge e){
		
		if (!H.edgeExists(v.getIndex(), w.getIndex())){ //no edge between nodes v and w
			
			//Find the minimum vertex cut between va and vb in modH
			MinVertexCut mvc = new MinVertexCut();
			mvc.findMinVertexCut(H, v.getIndex(), w.getIndex());	
			
			System.out.println("mvc: "+mvc.getCutSet().size()+"(L: "+mvc.getLset().size()+", R: "+mvc.getRset().size()+")");	
			
			//Return the M set for the computed separation into the (artificial) edge e
			e.setM(mvc.getCutSet());
			
			if ( (mvc.getCutSet().size()==2) && (mvc.getLset().size()>0) && (mvc.getRset().size()>0) ){ //a valid 2-separation
				
				//Compute the corresponding edge partition of the edges in D[]
				completeXY(mvc.getLset(),mvc.getRset(),D,sx,sy,false);
				System.out.println("sp: sx("+sx.size()+"), sy("+sy.size()+")");
				
				return true;
			}
		}
		
		return false;
	}
	
	//Determines if there is a k-separation of degree k between vertices v and w, 
	//		such that either the min vertex cut is an independent set in the original graph, 
	//		or the L and R sets of the separation each contain at least two nodes
	private boolean findKsepPairV(BranchTree H, BranchNode v, BranchNode w, int k, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy, BranchEdge D[], BranchEdge e){
		
		if (!H.edgeExists(v.getIndex(), w.getIndex())){ //no edge between nodes v and w
			
			//First, compute the min vertex cut between v and w
			MinVertexCut mvc = new MinVertexCut();
			mvc.findMinVertexCut(H, v.getIndex(), w.getIndex());
			System.out.println("mvc: "+mvc.getCutSet().size()+"(L: "+mvc.getLset().size()+", R: "+mvc.getRset().size()+")");
		
			//Get the M set from the graph
			e.setM(mvc.getCutSet());
			
			if (mvc.getCutSet().size()==k){ //found a k-separation
				
				if ( ((mvc.getLset().size()>=2)&&(mvc.getRset().size()>=2)) || (!isIndependentSetGraph(mvc.getCutSet())) ) { //both L and R have size at least two, or M is not an independent set in the original graph
				
					//Compute the corresponding edge partition of the edges in D[]
					completeXY(mvc.getLset(), mvc.getRset(), D, sx, sy, false);
					System.out.println("sp (nis): sx(" + sx.size() + "), sy(" + sy.size() + ")");
					
					//Get the M set from the graph
					e.setM(mvc.getCutSet());

					return true;
				}
				
				if ( (v.getNumEdges()<=k) && (w.getNumEdges()<=k) ) { //both v and w are of degree at most k
					
					//Find the sets of neighbors of nodes v and w
					LinkedHashSet<BranchNode> vN = H.getNeighborsForNode(v,false);
					LinkedHashSet<BranchNode> wN = H.getNeighborsForNode(w,false);
					
					Object vNa[] = vN.toArray();
					Object wNa[] = wN.toArray();
					
					//Contract v to each of its neighbors that are not also neighbors of w, and similarly for w,
					//		and check if any of the induced graphs generates a separation of size at most k
					for (int i=0; i<vNa.length; i++) { //for all neighbors of v
						
						if (!wN.contains((BranchNode)vNa[i])) { //current neighbor of v is not a neighbor of w
							
							for (int j=0; j<wNa.length; j++){ //for all neighbors of w
								
								if (!vN.contains((BranchNode)wNa[j])) { //current neighbor of w is not a neighbor of v
									
									if (!H.edgeExists(((BranchNode)vNa[i]).getIndex(),((BranchNode)wNa[j]).getIndex())) { //no edge between the current v and w neighbors
										
										int xyInd[] = new int[2];
										BranchTree modHvw = contractV(H,v,(BranchNode)vNa[i],vNa,w,(BranchNode)wNa[j],wNa,xyInd); //contract v and w
											
										//Compute the min vertex cut for the contracted graph between the new x and y nodes
										mvc = new MinVertexCut();
										mvc.findMinVertexCut(modHvw, xyInd[0], xyInd[1]);
										System.out.println("mvc (2c): "+mvc.getCutSet().size()+"(L: "+mvc.getLset().size()+", R: "+mvc.getRset().size()+")");
										
										if (mvc.getCutSet().size()<=k) { //found a k-separator with at least two nodes in each of the L and R sets of the original graph H
											
											//Add the contracted vertex names to the L (node v) and R (node w) sets
											mvc.addToLset(v.getv1());
											mvc.addToRset(w.getv1());
											
											//Compute the corresponding edge partition of the edges in D[]
											completeXY(mvc.getLset(), mvc.getRset(),D,sx,sy,false);
											System.out.println("sp: sx("+sx.size()+"), sy("+sy.size()+")");
											
											//Get the M set from the contracted graph
											e.setM(mvc.getCutSet());
											
											return true;
										}
									}
								}
							}
						}
					}
				}
				
				else if ( (v.getNumEdges()<=k) || (w.getNumEdges()<=k) ) { //one of v and w has degree at most k, and the other one has degree greater than k
					
					BranchNode sn = null; //the node of smaller degree (either v or w)
					BranchNode ln = null; //the node of larger degree (either v or w)
					if (v.getNumEdges()<=k) {
						sn = v;
						ln = w;
					}
					else {
						sn = w;
						ln = v;
					}
					
					LinkedHashSet<BranchNode> sN = H.getNeighborsForNode(sn,false);
					Object sNa[] = sN.toArray();
					
					//Contract sn to each of its neighbors,
					//		and check if any of the induced graphs generates a separation of size at most k
					for (int i=0; i<sNa.length; i++) { //for all neighbors of sn
						
						if (!H.edgeExists(((BranchNode)sNa[i]).getIndex(),ln.getIndex())) { //no edge between current neighbor of sn and node ln
							
							int xyInd[] = new int[2];
							BranchTree modHvw = contractV(H,sn,(BranchNode)sNa[i],sNa,xyInd,ln); //contract v
							
							//Compute the min vertex cut for the contracted graph between the new node x and the original node ln
							mvc = new MinVertexCut();
							mvc.findMinVertexCut(modHvw, xyInd[0], xyInd[1]);
							System.out.println("mvc (1c): "+mvc.getCutSet().size()+"(L: "+mvc.getLset().size()+", R: "+mvc.getRset().size()+")");
							
							if (mvc.getCutSet().size()<=k) { //found a k-separator with at least two nodes in each of the L and R sets of the original grpah H
								
								//Add the contracted vertex name (node sn) to the L set
								mvc.addToLset(sn.getv1());
								
								//Compute the corresponding edge partition of the edges in D[]
								completeXY(mvc.getLset(), mvc.getRset(),D,sx,sy,false);
								System.out.println("sp: sx("+sx.size()+"), sy("+sy.size()+")");
								
								//Get the M set from the contracted graph
								e.setM(mvc.getCutSet());
								
								return true;
							}
						}
					}
				}
				
				else { //both v and w have degrees greater than k
					
					//Compute the corresponding edge partition of the edges in D[]
					completeXY(mvc.getLset(), mvc.getRset(),D,sx,sy,false);
					System.out.println("sp: sx("+sx.size()+"), sy("+sy.size()+")");
					
					//Get the M set from the contracted graph
					e.setM(mvc.getCutSet());
					
					return true;
				}
			}
			
			else if (mvc.getCutSet().size()<k) {
				System.out.println("ERROR: looking for a "+k+"-separator, but a "+mvc.getCutSet().size()+"-separator found");
				System.exit(1);
			}
		}
		
		return false; //no 3-separation with the required properties found
	}
	
	//Checks if the set s of vertex names is an independent set in the original graph (i.e., if there are no edges in the graph between any two vertices in s)
	private boolean isIndependentSetGraph(LinkedHashSet<String> s){
		
		Object sA[] = s.toArray();
		for (int i=0; i<sA.length; i++){
			
			for (int j=i+1; j<sA.length; j++){
				
				for (int k=0; k<bt.getNumNodes(); k++) {
					
					BranchNode bn = bt.getNode(k);
					
					if (bn.getIsLeaf()) { //a leaf node in the partial branch decomposition, so this node corresponds to an edge in the original graph
						
						if ( (bn.getv1().equalsIgnoreCase((String)sA[i])) && (bn.getv2().equalsIgnoreCase((String)sA[j])) ) //edge exists between sA[i] and sA[j]
							return false;
						
						else if ( (bn.getv2().equalsIgnoreCase((String)sA[i])) && (bn.getv1().equalsIgnoreCase((String)sA[j])) ) //edge exists between sA[i] and sA[j]
							return false;
					}
				}
			}
		}
		return true; //no two vertices in the set are adjacent in H
	}
	
	//Checks if the set s of vertex names is an independent set in H (i.e., if there are no edges in H between any two vertices in s)
	private boolean isIndependentSet(BranchTree H, LinkedHashSet<String> s){
		
		Object sA[] = s.toArray();
		for (int i=0; i<sA.length; i++){
			
			int ind1 = getIndOfVertex(H,(String)sA[i]);
			
			for (int j=i+1; j<sA.length; j++){
				
				int ind2 = getIndOfVertex(H,(String)sA[j]);
				
				if (H.edgeExists(ind1, ind2)) //an edge exists between the nodes with indices ind1 and ind2, so the set is not independent
					return false;
			}
		}
		return true; //no two vertices in the set are adjacent in H
	}
	
	//Returns the index of the vertex with name vs into the node array of H
	private int getIndOfVertex(BranchTree H, String vs){
		
		for (int i=0; i<H.getNumNodes(); i++){
			if (H.getNode(i).getv1().equalsIgnoreCase(vs)) //found vertex (vertex name is stored in v1 of the current node)
				return H.getNode(i).getIndex();
		}
		
		System.out.println("ERROR: cannot find H node with vertex name "+vs);
		System.exit(1);
		return -1;
	}
	
	//Contracts node v to node x and node w to node y from graph h, and returns the modified graph;
	//NOTE: No changes to the input graph h are made here;
	//Node x in the modified graph is incident with the union of the edges for nodes x and v from the original graph, and node v is deleted;
	//Node y in the modified graph is incident with the union of the edges for nodes y and w from the original graph, and node w is deleted;
	//The nodes to which node v is adjacent in the original graph are given in vNa[] (a BranchNode array); similarly for w and wNa[]
	//The indices of nodes x and y in the modified graph are returned in xyInd[]
	private BranchTree contractV(BranchTree h, BranchNode v, BranchNode x, Object vNa[], BranchNode w, BranchNode y, Object wNa[], int xyInd[]){
		
		BranchTree modH = h.deepCopy(); //make a deep copy of the original graph
		
		//Get copies of the v, x, w, and y nodes for the modified graph  (before modifying the graph)
		BranchNode vMod = modH.getNode(v.getIndex());
		BranchNode xMod = modH.getNode(x.getIndex());
		BranchNode wMod = modH.getNode(w.getIndex());
		BranchNode yMod = modH.getNode(y.getIndex());
		
		//Add the new edges to nodes x and y in the modified graph
		unifyEdges(modH,vNa,xMod);
		unifyEdges(modH,wNa,yMod);		
		
		//Delete nodes v and w from the modified graph;
		//NOTE: this must happen after the edge copying
		modH.deleteNode(vMod.getIndex()); 
		modH.deleteNode(wMod.getIndex());
		
		//Return the indices of nodes x and y in the modified graph
		xyInd[0] = xMod.getIndex();
		xyInd[1] = yMod.getIndex();
		
		return modH;
	}
	
	//Contracts node v to node x  from graph h, and returns the modified graph;
	//NOTE: No changes to the input graph h are made here;
	//Node x in the modified graph is incident with the union of the edges for nodes x and v from the original graph, and node v is deleted;
	//The nodes to which node v is adjacent in the original graph are given in vNa[] (a BranchNode array);
	//The indices of nodes x and y in the modified graph are returned in xyInd[]
	private BranchTree contractV(BranchTree h, BranchNode v, BranchNode x, Object vNa[], int xyInd[], BranchNode ln){
		
		BranchTree modH = h.deepCopy(); //make a deep copy of the original graph
		
		//Get copies of the v and x nodes for the modified graph (before modifying the graph)
		BranchNode vMod = modH.getNode(v.getIndex());
		BranchNode xMod = modH.getNode(x.getIndex());
		BranchNode lnMod = modH.getNode(ln.getIndex());
		
		//Add the new edges to node x in the modified graph
		unifyEdges(modH,vNa,xMod);
		
		//Delete node v from the modified graph;
		//NOTE: this must happen after the edge copying
		modH.deleteNode(vMod.getIndex()); 
		
		//Return the indices of nodes x and y in the modified graph
		xyInd[0] = xMod.getIndex();
		xyInd[1] = lnMod.getIndex();
		
		return modH;
	}
	
	//Adds edges in modH between node x and all nodes in vNa[];
	//NOTE: this changes the input graph h
	private void unifyEdges(BranchTree h, Object vNa[], BranchNode x){
		
		for (int i=0; i<vNa.length; i++){ //add edges to all nodes in vNa to the set of edges for x in the graph h
			
			if (x.getIndex()!=((BranchNode)vNa[i]).getIndex()) //not the same node
				
				h.addEdge(x.getIndex(), ((BranchNode)vNa[i]).getIndex() );
		}
	}
	
	
	//Perform the eigenvector heuristic for finding edge partitions on node pn; the resulting partition is returned in the sets sx and sy
	private void doEigen(BranchNode pn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy){
		
		BranchEdge D[] = bt.getEdgesForNode(pn.getIndex()); //all the edges incident with node pn
		
		
		//First, construct the F matrix
		
		LinkedHashSet<String> S = new LinkedHashSet<String>(); //the union of the M sets for all edges in D[]
		for (int i=0; i<D.length; i++){
			S.addAll(D[i].getM());
		}		
		Object Sa[] = S.toArray();
		
		BranchEdge Nv[][] = new BranchEdge[Sa.length][D.length]; //the set of edges in D incident with each graph vertex in Sa
		for (int i=0; i<Sa.length; i++){
			int cnt = 0;
			for (int j=0; j<D.length; j++){
				if (D[j].getM().contains((String)Sa[i])){ //add edge D[a] to the list of adjacent edges in Nv[i]
					Nv[i][cnt] = D[j];
					cnt++;
				}
			}
			BranchEdge tmp[] = new BranchEdge[cnt];
			System.arraycopy(Nv[i], 0, tmp, 0, cnt);
			Nv[i] = tmp; //reduce array size to the actual number
		}
		
		double F[][] = new double[D.length][D.length]; //a |D|x|D| matrix
		for (int i=0; i<F.length; i++){
			for (int j=0; j<F[i].length; j++){
				F[i][j] = 0.0;
				if (i!=j){
					for (int k=0; k<Nv.length; k++){
						if ( isElement(Nv[k],D[i]) && isElement(Nv[k],D[j]) ) //both D[i] and D[j] are in Nv[k]
							F[i][j] += (-1.0)/(Nv[k].length-1);
					}
				}
				else { //diagonal entry
					for (int k=0; k<Nv.length; k++){
						if (isElement(Nv[k],D[i])) //D[i] is in Nv[k]
							F[i][j]++;
					}
				}
			}
		}
		
		
		//Find the eigenvalues and eigenvectors of F
		Matrix Fm = new Matrix(F);
		EigenvalueDecomposition ed = new EigenvalueDecomposition(Fm);
		double eVl[][] = ed.getD().getArray(); //the eigenvalues of F
		double eVc[][] = ed.getV().getArray(); //the eigenvectors of F
		
		//Find the second smallest eigenvalue of F
		double evlA[] = new double[eVl.length];
		for (int i=0; i<eVl.length; i++){
			evlA[i] = eVl[i][i];
		}
		evlA = sort(evlA); //sort in ascending order
		int evlind = -1;
		for (int i=0; i<eVl.length; i++){
			if (eVl[i][i]==evlA[1]) { //the second smallest eigenvalue
				evlind = i;
				break;
			}
		}
		
		
		//Find the eigenvector corresponding to the second smallest eigenvalue, and order D so that secondevc is in non-decreasing order
		double secondEvc[] = getColumnCopy(eVc, evlind);
		double secondEvcSorted[] = sort(secondEvc); //sort in non-decreasing order
		BranchEdge Dsorted[] = new BranchEdge[D.length];
		for (int i=0; i<secondEvcSorted.length; i++){
			for (int j=0; j<secondEvc.length; j++){
				if (secondEvc[j]==secondEvcSorted[i]){
					Dsorted[i] = D[j];
					secondEvc[j] = -Math.pow(10, 38);
					break;
				}
			}
		}
		
		//Add the first |D|/3 elements of the sorted D to the set sx, and the last |D|/3 elements to the set sy
		int numE = (int)Math.ceil(D.length/3.0);
		for (int i=0; i<numE; i++){
			sx.add(Dsorted[i]);
			sy.add(Dsorted[Dsorted.length-1-i]);
		}
		
		
		//Construct a new graph H with nodes the vertices in the union of the M sets for all edges in D[],
		//		and an edge between a pair of nodes if some M set for an edge in D[] contains both these nodes
		BranchTree H = constructGraphH(Sa,D);
		
		//Modify H by adding two vertices: va(source) and vb(sink); 
		//		add edges to va for all graph vertices incident with edges in sx;
		//		add edges to vb for all graph vertices incident with edges in sy;
		//		return the indices for va and vb of the modified graph into vi[]
		int vi[] = new int[2];
		H = modifyGraphH(H,sx,sy,vi,null,true);
		
		
		//Find the minimum vertex cut between va and vb in H
		MinVertexCut mvc = new MinVertexCut();
		mvc.findMinVertexCut(H.deepCopy(), vi[0], vi[1]);
		
		mvc.getCutSet().remove(""); //remove entries not corresponding to graph vertices
		mvc.getLset().remove("");
		mvc.getRset().remove("");
		
		System.out.println("mvc: "+mvc.getCutSet().size()+"(L: "+mvc.getLset().size()+", R: "+mvc.getRset().size()+")");
		
		//Get the computed Lset and Rset from mvc, and add all remaining edges of D to sx and sy (this gives the full edge partition)
		completeXY(mvc.getLset(),mvc.getRset(),D,sx,sy,false);
		
		System.out.println("sp: sx("+sx.size()+"), sy("+sy.size()+")");
	}
	
	//Determines if a[] contains e
	private boolean isElement(Object a[], Object e){
		
		for (int i=0; i<a.length; i++){
			if (a[i]==e)
				return true;
		}
		return false;
	}
	
	//Construct a graph H with nodes the vertices in the union of the M sets for all edges in D[] (stored in Sa[]),
	//		and an edge between a pair of nodes if some M set for an edge in D[] contains both these nodes
	private BranchTree constructGraphH(Object Sa[], BranchEdge D[]){
		
		BranchTree H = new BranchTree(); //construct a new graph
		for (int i=0; i<Sa.length; i++){ //add the graph vertices
			H.addNode(new BranchNode(true,(String)Sa[i],null)); //store the vertex name in v1 for the nodes
		}
		
		for (int i=0; i<H.getNumNodes(); i++){ //add the edges
			
			String v1 = H.getNode(i).getv1(); //the vertex name is stored in v1 for the nodes
			
			for (int j=i+1; j<H.getNumNodes(); j++){
				
				String v2 = H.getNode(j).getv1(); //the vertex name is stored in v1 for the nodes
				
				for (int k=0; k<D.length; k++){
					
					if ( D[k].getM().contains(v1) && D[k].getM().contains(v2) ){
						
						H.addEdge(i,j); //add the edge
						break;
					}
				}
			}
		}
		
		return H;
	}
	
	//Modify the graph H by adding to vertices: va(source) and vb(sink); 
	//		add edges to va for all graph vertices incident with edges in the set sx;
	//		add edges to vb for all graph vertices incident with edges in the set sy;
	//		if (sy==null), then there is a single sink specified by the singleSink node;
	//		if (eigen==true), then the set sy consists of edges in the transformed graph; otherwise, the set sy consists of edges in the original graph H
	private BranchTree modifyGraphH(BranchTree H, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy, int vi[], BranchNode singleSink, boolean eigen){
		
		BranchNode va = new BranchNode(true,"",null); //all vertices in the M sets in sx are identified to va
		H.addNode(va);		
		if (eigen)
			modifyGraphHhelperEigen(H,va,sx);
		else
			modifyGraphHhelper23sep(H,va,sx);
		
		BranchNode vb = null;
		if (sy!=null) { //multiple sinks are present
			vb = new BranchNode(true,"",null); //all vertices in the M sets in sy are identified to vb
			H.addNode(vb);
			if (eigen)
				modifyGraphHhelperEigen(H,vb,sy);
			else
				modifyGraphHhelper23sep(H,vb,sy);
		}
		
		vi[0] = va.getIndex(); //return the indices for va and vb in the modified graph
		if (sy!=null)
			vi[1] = vb.getIndex();
		else
			vi[1] = singleSink.getIndex();
		
		return H;
	}
	
	//Helper for modifyGraphH(); adds edges to node vn for all vertices in the M sets of the edges in set s
	private void modifyGraphHhelperEigen(BranchTree H, BranchNode vn, LinkedHashSet<BranchEdge> s){
		
		Object sA[] = s.toArray();
		LinkedHashSet<String> nodes = new LinkedHashSet<String>();
		for (int i=0; i<sA.length; i++){
			Object curM[] = ((BranchEdge)sA[i]).getM().toArray();
			for (int j=0; j<curM.length; j++){
				nodes.add((String)curM[j]);
			}
		}
		modifyGraphHaddEdges(H,vn,nodes);
	}
	
	//Helper for modifyGraphH(); adds edges to node vn for all vertices incident with the edges in set s
	private void modifyGraphHhelper23sep(BranchTree H, BranchNode vn, LinkedHashSet<BranchEdge> s){
		
		Object sA[] = s.toArray();
		LinkedHashSet<String> nodes = new LinkedHashSet<String>();
		for (int i=0; i<sA.length; i++){
			nodes.add(((BranchEdge)sA[i]).getn1().getv1()); //the vertex name is stored in v1 of the nodes
			nodes.add(((BranchEdge)sA[i]).getn2().getv1());
		}
		modifyGraphHaddEdges(H,vn,nodes);
	}
	
	//Adds an edge between node vn and each node in the set nodes in graph H
	private void modifyGraphHaddEdges(BranchTree H, BranchNode vn, LinkedHashSet<String> nodes){
		
		Object nodesA[] = nodes.toArray();
		for (int i=0; i<nodesA.length; i++){ //add edges between va and the nodes in H for each vertex in vaNodes
			for (int j=0; j<H.getNumNodes(); j++){
				if ( ((String)nodesA[i]).equalsIgnoreCase(H.getNode(j).getv1()) ){
					H.addEdge(vn.getIndex(),H.getNode(j).getIndex()); //add edge
					break;
				}
			}
		}
	}
	
	//Compute the complete partition of the set of edges in D[] into the two edge sets sx and sy;
	//For each edge e in D[], 
	//		if the intersection of the set M of e with Rset is empty and the the intersection of the set M of e with Lset is empty, 
	//				add e to X (if inX==true) or to Y (otherwise);
	//		otherwise, if the intersection of the set M of e with Rset is empty, then add e to X;
	//		otherwise, if the intersection of the set M of e with Lset is empty, then add e to Y;
	//		otherwise, output an error message
	private void completeXY(LinkedHashSet<String> Lset, LinkedHashSet<String> Rset, BranchEdge D[],
			LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy, boolean inX){
		
		for (int i=0; i<D.length; i++){
			
			if ( ! (sx.contains(D[i]) || sy.contains(D[i]) ) ){ //edge not already assigned to one of the two sets
			
				LinkedHashSet<String> eM = D[i].getM();
				
				LinkedHashSet<String> intersectL = new LinkedHashSet<String>(Lset);
				intersectL.retainAll(eM);
				
				LinkedHashSet<String> intersectR = new LinkedHashSet<String>(Rset);
				intersectR.retainAll(eM);
				
				if ( intersectR.isEmpty() && intersectL.isEmpty() ){
					if (inX)
						sx.add(D[i]);
					else
						sy.add(D[i]);
				}
				
				else if (intersectR.isEmpty())
					sx.add(D[i]);
				
				else if (intersectL.isEmpty())
					sy.add(D[i]);
				
				else {
					System.out.println("ERROR: current edge cannot be assigned to a partition");
					System.exit(1);
				}
			}
		}		
	}
	
	//Performs the split at the internal node bn
	private void performSplit(BranchNode bn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy){
		
		//Split bn into two nodes (bx and by), such that all edges in sx are adjacent to bx and all edges
		//		in sy are adjacent to by, and there is an edge between bx and by
		
		BranchNode bx = new BranchNode(false,null,null);
		bt.addNode(bx);
		
		BranchNode by = new BranchNode(false,null,null);
		bt.addNode(by);		
		
		BranchEdge bnEdges[] = bt.getEdgesForNode(bn.getIndex());
		
		for (int i=0; i<bnEdges.length; i++){
			
			BranchNode otherNode = null;			
			if (bnEdges[i].getn1().getIndex()==bn.getIndex())
				otherNode = bnEdges[i].getn2();
			else
				otherNode = bnEdges[i].getn1();
			
			int eInd = -1;			
			if (sx.contains(bnEdges[i])) //current edge is in sx				
				eInd = bt.addEdge(bx.getIndex(), otherNode.getIndex());
			
			else //current edge is in sy
				eInd = bt.addEdge(by.getIndex(), otherNode.getIndex());
			
			BranchEdge ce = bt.getEdge(eInd); //copy the M set to the new edge
			ce.setM(bnEdges[i].getM());
			
			ce.setChecked(otherNode, bnEdges[i].isChecked(otherNode)); //copy the is-checked flag for otherNode from the original edge
			if (sx.contains(bnEdges[i])) //original edge is in sx
				ce.setChecked(bx, bnEdges[i].isChecked(bn)); //copy the is-checked flag for node bn from the original edge
			else //original edge is in sy
				ce.setChecked(by, bnEdges[i].isChecked(bn));
		}
		
		bt.deleteNode(bn.getIndex()); //delete the original node
		
		int bxbyInd = bt.addEdge(bx.getIndex(),by.getIndex()); //add an edge between the new nodes bx and by
		BranchEdge bxby = bt.getEdge(bxbyInd);
		
		compMLRsets(bxby,sx,sy); //compute the M, L, and R sets for the edge (bx-by)
		
		System.out.println("done.. Width of new edge: "+bxby.getM().size());
	}
	
	//Compute the M, L, and R sets for edge e, such that all edges in the set sx are incident with one of the endpoints of e,
	//		while all edges in the set sy are incident with the other endpoint of edge e;
	//Assumes that the M, L, and R sets for edges in the sets sx and sy have already been computed
	private void compMLRsets(BranchEdge e, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy){
		
		Object sxa[] = sx.toArray();
		LinkedHashSet<String> uMsx = new LinkedHashSet<String>(); //the union of the graph vertices in the M sets of all edges in sx
		for (int i=0; i<sxa.length; i++)
			uMsx.addAll(((BranchEdge)sxa[i]).getM());
		
		Object sya[] = sy.toArray();
		LinkedHashSet<String> uMsy = new LinkedHashSet<String>(); //the union of the graph vertices in the M sets of all edges in sy
		for (int i=0; i<sya.length; i++)
			uMsy.addAll(((BranchEdge)sya[i]).getM());
		
		LinkedHashSet<String> xyM = new LinkedHashSet<String>(uMsx);
		xyM.retainAll(uMsy); //the intersection of the two set unions (uMsx and uMsy)
		e.setM(xyM); //the M set of edge (bx-by)
	}
	
	//Reads in the residue interaction graph from graphFile;
	//The set gvDup returns all graph vertices that are part of more than one node (this is used for initializing the star graph)
	private void readGraphFile(String graphFile, LinkedHashSet<String> gvDup){
		
		bt = new BranchTree();
		gv = new GraphVertices();
		
		BufferedReader bufread = null;
		try {
			File file = new File(graphFile);		
			FileReader fr = new FileReader(file);  
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			System.out.println("ERROR: Residue interaction graph file not found");
			System.exit(1);
		}
		
		String str = null;
		boolean done = false;
		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
				System.out.println("ERROR: An error occurred while reading input");
				System.exit(1);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			
			else if ( (!getToken(str,1).equalsIgnoreCase("PIG:0")) && (!str.equalsIgnoreCase("0 0")) ) { //not the first/last lines in file
				
				String v1 = getToken(str,1);
				String v2 = getToken(str,2);
				
				bt.addNode(new BranchNode(true,v1,v2)); //add the new node to bt
				
				if (!gv.addV(v1)) //add to the graphVertices
					gvDup.add(v1); //vertex already found in graph, so update gvDup
				
				if (!gv.addV(v2))
					gvDup.add(v2);
			}
		}
		
		gv.setCompleteGraph(); //all graph vertices have been added
		
		try { bufread.close(); } catch(Exception e){} //done, so close file for reading
	}
	
	//Output the computed branch decomposition
	private void outputBranchDecomposition(String fout){
		
		//Output branch decomposition for viewing
		PrintStream bd = setupOutputFile(fout+"_v.txt");
		bd.println("PIG:0 "+fout);		
		for (int i=0; i<bt.getNumEdges(); i++){			
			BranchEdge be = bt.getEdge(i);		
			bd.println(be.getn1().getIndex()+" "+be.getn2().getIndex());
		}
		bd.println("0 0");		
		bd.flush();
		bd.close();
		
		
		//Output branch decomposition to be read as input to the protein design algorithm
		bd = setupOutputFile(fout);
		
		//First, output node information
		bd.println("NODES: "+bt.getNumNodes());
		bd.println("# NODE LINE FORMAT:");
		bd.println("# node_index num_edges is_leaf graph_v1 graph_v2");
		
		for (int i=0; i<bt.getNumNodes(); i++){
			
			BranchNode bn = bt.getNode(i);
			
			bd.print(bn.getIndex()+" "+bn.getNumEdges()+" "+bn.getIsLeaf());
			if (bn.getIsLeaf())
				bd.print(" "+bn.getv1()+" "+bn.getv2());
			bd.println();
		}
		bd.println();
		
		
		//Output edge information
		bd.println("EDGES: "+bt.getNumEdges());
		bd.println("# EDGE LINES FORMAT:");
		bd.println("# edge_index node1_index node2_index width");
		bd.println("# list_graph_vertices_in_M");
		
		int branchWidth = 0;
		for (int i=0; i<bt.getNumEdges(); i++){
			
			BranchEdge be = bt.getEdge(i);
			BranchNode bn1 = be.getn1();
			BranchNode bn2 = be.getn2();
			
			Object curM[] = be.getM().toArray();
			
			branchWidth = Math.max(branchWidth,curM.length);
			
			bd.println(i+" "+bn1.getIndex()+" "+bn2.getIndex()+" "+curM.length);
			
			if (curM.length<=0){
				bd.print("-1");
			}
			else {
				for (int j=0; j<curM.length; j++){
					bd.print(curM[j]+" ");
				}
			}
			bd.println();
		}
		bd.println();
		bd.println("# Branchwidth: "+branchWidth);
		bd.flush();
		bd.close();
	}
	
	//Setup the file with name filename for output
	private PrintStream setupOutputFile(String fileName){
		PrintStream logPS = null; //the output file for conf info
		try {			
			FileOutputStream fileOutputStream = new FileOutputStream(fileName);
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			logPS = new PrintStream( bufferedOutputStream );
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
		}
		return logPS;
	}
	
	//Sorts a[] in non-decreasing order
	private double [] sort(double a[]){
		
		double aSorted[] = new double[a.length];
		System.arraycopy(a, 0, aSorted, 0, aSorted.length);
		Arrays.sort(aSorted);
		return aSorted;
	}
	
	//Returns a copy of the elements in column cInd of the square array a[][]
	private double [] getColumnCopy(double a[][], int cInd){
		
		double c[] = new double[a.length];
		for (int i=0; i<c.length; i++){
			c[i] = a[i][cInd];
		}
		return c;
	}
	
	// This function returns the xth token in string s
	private String getToken(String s, int x) {
	
		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");
		
		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
			  st.nextToken();
			else {
				// System.out.println("ERROR: Unable to access parameter " + x);
				return(new String(""));
			}
		}
		
		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken
	
	//Generates combinations for a set of n elements, such that each element is either in the '0' or the '1' set,
	//		and such that both the '0' and '1' sets are non-empty;
	//		symmetric combinations (e.g., 0011 and 1100) are not counted as distinct combinations;
	//Returns the set of combinations in comb[][]
	private int [][] genComb(int n){
		
		int numComb = (int)(Math.pow(2, n)-2)/2;
		int comb[][] = new int[numComb][n];
		
		int curNum[] = new int[1];
		curNum[0] = 0;
		
		int curComb[] = new int[n];
		
		genCombHelper(n,comb,curNum,curComb,0);
		
		return comb;
	}
	
	//Helper for genComb(); generates all desired combinations and returns them in comb[][]
	private void genCombHelper(int n, int comb[][], int curNum[], int curComb[], int depth){
		
		if (depth>=n){ //full combination, so assign
			
			for (int i=0; i<comb[curNum[0]].length; i++) {
				comb[curNum[0]][i] = curComb[i];
			}
			curNum[0]++;
		}
		else { //assign a set for the current depth
			
			if ( !( (depth==(n-1)) && (curNum[0]==0) ) ){ //do not generate the set of all 0's			
				curComb[depth] = 0; //first, generate all combinations with the current element in the '0' set
				genCombHelper(n,comb,curNum,curComb,depth+1);
			}
			
			if (depth>0){
				curComb[depth] = 1; //then, generate all combinations with the current element in the '1' set
				genCombHelper(n,comb,curNum,curComb,depth+1);
			}
		}
	}
	
	//Generate all pairwise combinations for the elements in a[]
	private int [][] genPairComb(int a[]) {
		
		int numComb = a.length*(a.length-1) / 2;
		int comb[][] = new int[numComb][2];
		
		int curComb = 0;
		for (int i=0; i<a.length; i++){
			for (int j=i+1; j<a.length; j++){
				comb[curComb][0] = a[i];
				comb[curComb][1] = a[j];
				curComb++;
			}
		}
		
		return comb;
	}
	
	//Functions for calculating CPU Time
	public long CPUTime(){
		ThreadMXBean thread = ManagementFactory.getThreadMXBean();
		if(thread.isCurrentThreadCpuTimeSupported())
			return thread.getCurrentThreadCpuTime();
		else
			return 0L;
	}
	public long UserTime(){
		ThreadMXBean thread = ManagementFactory.getThreadMXBean();
		if(thread.isCurrentThreadCpuTimeSupported())
			return thread.getCurrentThreadUserTime();
		else
			return 0L;
	}
	public long SystemTime(){

		ThreadMXBean thread = ManagementFactory.getThreadMXBean();
		if(thread.isCurrentThreadCpuTimeSupported())
			return thread.getCurrentThreadCpuTime() - thread.getCurrentThreadUserTime();
		else
			return 0L;
	}
}
