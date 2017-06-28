package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.sparse.Subproblem;
import edu.duke.cs.osprey.sparse.TreeEdge;
import edu.duke.cs.osprey.tools.ResidueIndexMap;

public class TestSubproblem extends Subproblem {

	public TestSubproblem (RCs superSpace, TreeEdge sparseTree, ResidueIndexMap resMap, RCTuple initialConf) {
		super(superSpace, sparseTree, resMap, initialConf);
	}
	
	@Override
	public void preprocess()
	{
		int[] currentConf = new int[localConfSpace.getNumPos()];
		recursivelyProcessTuples(0,currentConf);
	}

}
