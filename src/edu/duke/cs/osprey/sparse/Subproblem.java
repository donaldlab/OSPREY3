package edu.duke.cs.osprey.sparse;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Map;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.tools.ResidueIndexMap;

public class Subproblem {
	private Set<Integer> LSet;
	private Set<Integer> lambdaSet;
	private Set<Integer> MSet;
	private Set<Integer> MULambdaSet;
	List<ConformationProcessor> processors = new ArrayList<>();
	public Subproblem leftSubproblem;
	public Subproblem rightSubproblem;
	private ResidueIndexMap residueIndexMap;
	protected RCs localConfSpace;
	private Map<Integer, Map<Integer, Integer>> PDBRCToSubspaceRCMap;
	private Map<Integer, Map<Integer, Integer>> SubspaceRCToPDBRCMap;
	
	public Subproblem (RCs superSpace, TreeEdge sparseTree, ResidueIndexMap resMap) {
		this(superSpace, sparseTree, resMap, null);
	}
	
	public Subproblem (RCs superSpace, TreeEdge sparseTree, ResidueIndexMap resMap, RCTuple initialConf) {
		localConfSpace = superSpace;
		residueIndexMap = resMap;
		
		TreeEdge curEdge = sparseTree;
		lambdaSet = curEdge.getLambda();
		MSet = curEdge.getM();
		MULambdaSet = new HashSet<Integer>();
		MULambdaSet.addAll(lambdaSet);
		MULambdaSet.addAll(MSet);
		initSubSpaceRCMap();
		if(initialConf != null)
			localConfSpace = superSpace.returnSubspace(initialConf);

		generateSubproblems(sparseTree);
		
	}

	private void initSubSpaceRCMap () {
		PDBRCToSubspaceRCMap = new HashMap<>();
		SubspaceRCToPDBRCMap = new HashMap<>();
		for(int PDBIndex : MULambdaSet)
		{
			int RCSpaceIndex = residueIndexMap.PDBIndexToDesignIndex(PDBIndex);
			if(localConfSpace == null)
			{
				System.out.println("Null.");
			}
			for(int RCSpaceRCIndex = 0; RCSpaceRCIndex < localConfSpace.getNum(RCSpaceIndex); RCSpaceRCIndex++)
			{
				int PDBSpaceRCIndex = localConfSpace.get(RCSpaceIndex, RCSpaceRCIndex);
				mapPDBRCToRCSpace(PDBIndex, PDBSpaceRCIndex, RCSpaceIndex, RCSpaceRCIndex);
			}
		}
	}

	private void mapPDBRCToRCSpace (int PDBIndex, int PDBSpaceRCIndex, int RCSpaceIndex, int RCSpaceRCIndex) {
		if(!PDBRCToSubspaceRCMap.containsKey(PDBIndex))
			PDBRCToSubspaceRCMap.put(PDBIndex, new HashMap<>());
		PDBRCToSubspaceRCMap.get(PDBIndex).put(PDBSpaceRCIndex, RCSpaceRCIndex);		
		if(!SubspaceRCToPDBRCMap.containsKey(RCSpaceIndex))
			SubspaceRCToPDBRCMap.put(RCSpaceIndex, new HashMap<>());
		SubspaceRCToPDBRCMap.get(RCSpaceIndex).put(RCSpaceRCIndex, PDBSpaceRCIndex);
		
	}

	private int mapSubproblemConfToIndex(int[] localConf)
	{
		int subproblemIndex = 0;
		int multiplier = 1;
		for(int i = 0; i < localConf.length; i++)
		{
			int PDBIndex = residueIndexMap.designIndexToPDBIndex(i);
			int subspaceRCIndex = getRCSpaceRCIndex(PDBIndex, localConf[i]);
			subproblemIndex += subspaceRCIndex*multiplier;
			multiplier *= localConfSpace.getNum(i);
		}
		
		return subproblemIndex;
	}
	
	private int getRCSpaceRCIndex(int PDBIndex, int PDBSpaceRCIndex)
	{
		return PDBRCToSubspaceRCMap.get(PDBIndex).get(PDBSpaceRCIndex);
	}
	
	private void generateSubproblems (TreeEdge sparseTree) {


		TreeEdge curEdge = sparseTree;
		//ordered array for recursion
		if(curEdge.leftChild != null)
		{
			leftSubproblem = new Subproblem(localConfSpace, sparseTree.leftChild.getCofEdge(), residueIndexMap);
		}
		if(curEdge.rightChild != null)
		{	
			rightSubproblem = new Subproblem(localConfSpace, sparseTree.rightChild.getCofEdge(), residueIndexMap);
		}
	}



	public void addConformationProcessor(ConformationProcessor processor)
	{
		if(leftSubproblem != null)
			leftSubproblem.addConformationProcessor(processor);
		if(rightSubproblem != null)
			rightSubproblem.addConformationProcessor(processor);
		processors.add(processor);
	}

	public BigInteger getTotalConformations()
	{
		return localConfSpace.unprunedConfsFromRCs();
	}
	
	public BigInteger getTotalLocalConformations()
	{
		BigInteger numConformations = BigInteger.ONE;
		for(int PDBIndex : MULambdaSet)
		{
			int designIndex = residueIndexMap.PDBIndexToDesignIndex(PDBIndex);
			numConformations = numConformations.multiply(BigInteger.valueOf(localConfSpace.get(designIndex).length));
		}
		return numConformations;
	}
	
	public BigInteger getSubtreeTESS()
	{
		BigInteger numConformations = getTotalLocalConformations();
		numConformations = numConformations.add(leftSubproblem.getSubtreeTESS());
		numConformations = numConformations.add(rightSubproblem.getSubtreeTESS());
		return numConformations;
	}

	public void preprocess () {
		if(leftSubproblem != null)
			leftSubproblem.preprocess();
		if(rightSubproblem != null)
			rightSubproblem.preprocess();
		int[] currentConf = new int[localConfSpace.getNumPos()];
		recursivelyProcessTuples(0,currentConf);
	}



	protected void recursivelyProcessTuples (int position, int[] currentConf) {
		if(position >= localConfSpace.getNumPos())
		{
			System.out.println("Process conformation:"+printConf(currentConf));
			RCTuple confTuple = new RCTuple(currentConf);
			for(ConformationProcessor proc : processors)
			{
				proc.processConformation(confTuple);
			}
			return;
		}
		
		int PDBIndex = residueIndexMap.designIndexToPDBIndex(position);
		if(!MULambdaSet.contains(PDBIndex))
		{
			recursivelyProcessTuples(position+1, currentConf);
			return;
		}
		
		for(int i = 0; i < localConfSpace.getNum(position); i++)
		{
			currentConf[position] = localConfSpace.get(position, i);
			recursivelyProcessTuples(position+1, currentConf);
		}
	}

	private String printSubspaceConf (int[] currentConf) {
		String output = "(";
		for(int i = 0; i < currentConf.length-1; i++)
		{
			output+=i+":"+currentConf[i]+", ";
		}
		output = output+(currentConf.length-1)+":"+currentConf[currentConf.length-1]+")";
		return output;
	}

	private String printConf (int[] currentConf) {
		String output = "(";
		for(int i = 0; i < currentConf.length-1; i++)
		{
			int PDBIndex = residueIndexMap.designIndexToPDBIndex(i);
			if(MULambdaSet.contains(PDBIndex))
				output+=PDBIndex+":"+currentConf[i]+", ";
		}
		int finalPDBIndex = residueIndexMap.designIndexToPDBIndex(currentConf.length-1);
		if(MULambdaSet.contains(finalPDBIndex))
			output += (finalPDBIndex)+":"+currentConf[currentConf.length-1];
		output += ")";
		return output;
	}

}
