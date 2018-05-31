package edu.duke.cs.osprey.sparse;

import edu.duke.cs.osprey.confspace.RCTuple;

public class SparseKStarScoreEvaluator{

	private ConformationProcessor KStarScoreComputer;
	private BranchDecomposedProblem searchProblem;

	public SparseKStarScoreEvaluator(BranchDecomposedProblem problem)
	{
		searchProblem = problem;
		initializeSubevaluators(problem.getRoot());
	}

	private void initializeSubevaluators (Subproblem problem) {
		ConformationProcessor scoreComputer = new PartialScoreComputer();
		problem.addConformationProcessor(scoreComputer);
		if(problem.leftSubproblem != null)
			initializeSubevaluators(problem.leftSubproblem);
		if(problem.rightSubproblem != null)
			initializeSubevaluators(problem.rightSubproblem);
	}
	
	public double computeSparseKStarScore(RCTuple sequence)
	{
		return 0;
	}

}
