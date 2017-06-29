package edu.duke.cs.osprey.sparse;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.LinkedHashSet;
import java.util.StringTokenizer;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.structure.Molecule;
public class BranchTree {
	
	ResidueInteractionGraph G;
	private edu.duke.cs.osprey.sparse.TreeNode root;
	
	public BranchTree(String fName, SearchProblem problem){
		
		int numV = problem.confSpace.numPos;
		
		G = new ResidueInteractionGraph();
		
		readBranchDecomposition(fName, problem.confSpace.m);
		traverseTree();
	}

	//Read the branch decomposition from file fName
	private void readBranchDecomposition(String fName, Molecule m){
		
		BufferedReader bufread = null;
		try {
			File file = new File(fName);		
			FileReader fr = new FileReader(file);  
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			System.out.println("ERROR: Residue interaction graph file not found");
			System.exit(1);
		}
		
		//Read the node information
		String str = readLine(bufread,fName);
		TreeNode tn[] = new TreeNode[new Integer(getToken(str,2)).intValue()];
		readLine(bufread,fName); //skip two comment lines
		readLine(bufread,fName);
		for (int i=0; i<tn.length; i++){ //read in the information for each node
			
			str = readLine(bufread,fName); //read node line
			
			int tnName = new Integer(getToken(str,1)).intValue();
			int tnEdges = new Integer(getToken(str,2)).intValue();
			boolean tnIsLeaf = new Boolean(getToken(str,3)).booleanValue();
			
			if (tnIsLeaf) { //leaf node
				int tnV1 = m.getResByPDBResNumber(getToken(str,4)).getPDBIndex(); //transform to PDB numbering
				int tnV2 = m.getResByPDBResNumber(getToken(str,5)).getPDBIndex();
				
				tn[i] = new TreeNode(tnName,tnIsLeaf,tnV1,tnV2); //add node
				
				//Add the two graph vertices and the corresponding edge to the residue interaction graph G
				G.addVertex(tnV1);
				G.addVertex(tnV2);
				G.addEdge(tnV1,tnV2);
			}
			
			else //internal node
				tn[i] = new TreeNode(tnName,tnIsLeaf,-1,-1);
		}
		
		//Read the edge information
		readLine(bufread,fName); //skip blank line
		str = readLine(bufread,fName);
		TreeEdge te[] = new TreeEdge[new Integer(getToken(str,2)).intValue()];
		readLine(bufread,fName); //skip three comment lines
		readLine(bufread,fName);
		readLine(bufread,fName);
		for (int i=0; i<te.length; i++){
			
			str = readLine(bufread,fName); //read edge line
			//int teName = new Integer(getToken(str,1)).intValue();
			int teNode1 = new Integer(getToken(str,2)).intValue();
			int teNode2 = new Integer(getToken(str,3)).intValue();
			int teWidth = new Integer(getToken(str,4)).intValue();
			
			str = readLine(bufread,fName); //read the M set
			LinkedHashSet<Integer> teM = new LinkedHashSet<Integer>();
			for (int j=0; j<teWidth; j++){
				teM.add(m.getResByPDBResNumber(getToken(str,j+1)).getPDBIndex());
			}
			
			te[i] = new TreeEdge(teNode1, teNode2, teM, false);
		}

		//Transform the branch decomposition into a rooted tree
		transformRootedTree(tn, te);

		try { bufread.close(); } catch(Exception e){} //done, so close file for reading
	}
	
	//Called by readBranchDecomposition();
	//Transforms the branch decomposition (nodes in tn[], edges in te[]) into a rooted tree:
	//		Select an arbitrary edge and add node s between the vertices incident with this edge; add the root and connect it to node s
	private void transformRootedTree(TreeNode tn[], TreeEdge te[]){
		
		//First, add the root and a node s, and the edge between root and s
		root = new TreeNode(-3,false,-1,-1);
		TreeNode s = new TreeNode(-2,false,-1,-1);
		TreeEdge rs = new TreeEdge(root.getName(),s.getName(),(new LinkedHashSet<Integer>()), true);
		setPCvar(root,s,rs,true); //set s to be the left (and only) child of the root
		
		
		//Next, select the first edge from te[], and use it to root the tree; delete the edge from te[]
		final int useEdgeInd = 0;
		TreeEdge se = te[useEdgeInd].deepCopy();
		TreeEdge tte[] = new TreeEdge[te.length-1];
		int cnt = 0;
		for (int i=0; i<te.length; i++) {
			if (i!=useEdgeInd) {
				tte[cnt] = te[i];
				cnt++;
			}
		}
		te = tte;
		
		int lcInd = getInd(tn,se.getNodeName1()); //set the left child of s 
		setPCvar(s,tn[lcInd],new TreeEdge(s.getName(),tn[lcInd].getName(),se.getM(), false),true);
		
		int rcInd = getInd(tn,se.getNodeName2()); //set the right child of s
		setPCvar(s,tn[rcInd],new TreeEdge(s.getName(),tn[rcInd].getName(),se.getM(), false),false);
		
		//Finally, build the tree by following the edges in te[]
		buildRootedTreeHelper(tn[lcInd],te,tn);
		buildRootedTreeHelper(tn[rcInd],te,tn);
	}
	
	//Find the edges incident with node r, update the corresponding variables, and follow the edges to build the entire tree
	private void buildRootedTreeHelper(TreeNode r, TreeEdge te[], TreeNode tn[]){
		
		if (r.isLeaf()) { //current node is leaf node, so only one incident edge (leading to this node)
			return;
		}
		
		else { //current node is internal, so exactly three incident edges (one leading to this node, and two starting from this node)
			
			//get the edges starting from node r, and the corresponding left and right children of r
			int cNames[] = new int[2];
			TreeEdge ce[] = new TreeEdge[2];
			getChildrenOfNode(r,te,cNames,ce);
			
			int lcInd = getInd(tn,cNames[0]); //set the left child of r
			setPCvar(r,tn[lcInd],ce[0],true);
			
			int rcInd = getInd(tn,cNames[1]); //set the right child of r
			setPCvar(r,tn[rcInd],ce[1],false);
			
			//recursively build the tree
			buildRootedTreeHelper(tn[lcInd],te,tn);
			buildRootedTreeHelper(tn[rcInd],te,tn);
		}
	}
	
	//Find the two children edges of the internal node r and the corresponding node names of the two children;
	//		do not call this when r is a leaf node
	private void getChildrenOfNode(TreeNode r, TreeEdge te[], int cNames[], TreeEdge ce[]){
		
		int cn = 0;
		for (int i=0; i<te.length; i++){
			
			if ( te[i].getNodeName1()==r.getName() ) { //found an edge incident with node r
				
				if ((te[i].getNodeName2()!=r.getp().getName())) { //the other node incident with this edge is not the parent of node r
					
					cNames[cn] = te[i].getNodeName2();
					ce[cn] = te[i];
					cn++;
				}
			}
			
			else if ( te[i].getNodeName2()==r.getName() ) { //found an edge incident with node r
				
				if ((te[i].getNodeName1()!=r.getp().getName())) { //the other node incident with this edge is not the parent of node r
					
					cNames[cn] = te[i].getNodeName1();
					ce[cn] = te[i];
					cn++;
				}
			}
			
			if (cn==cNames.length) //found both children
				break;
		}
		
		if (cn!=2) {
			System.out.println("ERROR: could not find all nodes adjacent with node "+r.getName());
			System.exit(1);
		}
	}
	
	//Sets the necessary variables for the parent node p, child node c (left if isL==true, right otherwise), and the edge pc between p and c
	private void setPCvar(TreeNode p, TreeNode c, TreeEdge pc, boolean isL){
		
		pc.setP(p);
		pc.setC(c);
		
		c.setP(p);
		c.setCofEdge(pc);
		
		if (isL)
			p.setLc(c);		
		else
			p.setRc(c);
	}
	
	//Returns the index of the element of a[] with name s; returns -1 if s is not found in a[]
	private int getInd(TreeNode a[], int s){
		for (int i=0; i<a.length; i++){
			if (a[i].getName()==s)
				return i;
		}
		return -1;
	}
	
	//Reads the next line from bufread
	private String readLine(BufferedReader bufread, String fName){
		try {
			return bufread.readLine();
		}
		catch ( Exception e ){
			System.out.println("ERROR: An error occurred while reading input from file "+fName);
			System.exit(1);
		}
		return null;
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


	public TreeNode getRoot() {
		return root.getlc();
	}
	
	public TreeEdge getRootEdge() {
		return root.getlc().getCofEdge();
	}


	public ResidueInteractionGraph getGraph() {
		return G;
	}

	//Traverses the current tree starting at the root
	public void traverseTree(){
		
		traverseTreeHelper(root);
	}

	//Performs post-order traversal of the tree
	private void traverseTreeHelper(TreeNode n){
		
		TreeNode lc = n.getlc();
		if (lc!=null) //traverse left subtree first
			traverseTreeHelper(lc);
		
		TreeNode rc = n.getrc();
		if (rc!=null) //then traverse right subtree
			traverseTreeHelper(rc);
		
		if (n==root) //done
			return;
		
		else { //not at the root, so do the computation for the edge for which node n is the child
			TreeEdge curEdge = n.getCofEdge();
			if (curEdge.getIsRootEdge()){ //this is the root edge, so finish computation and output results
				
				System.out.println("Starting A matrix computation for root edge..");
				curEdge.compLlambda(); //compute the L and lambda sets for the root edge
								
			}
			
			else {
				curEdge.compLlambda(); //compute the L and lambda sets for the current edge
				if (curEdge.getIsLambdaEdge()){ //lambda edge, so perform the A matrix computation
					System.out.println("lambda: "+curEdge.getLambda().size()+", L: "+curEdge.getL().size()+", M: "+curEdge.getM().size());
				}
			}
		}
	}
	
}
