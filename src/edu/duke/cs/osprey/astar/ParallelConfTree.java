package edu.duke.cs.osprey.astar;

import java.util.ArrayList;

import edu.duke.cs.osprey.astar.parallel.WorkCrew;
import edu.duke.cs.osprey.astar.parallel.Worker;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class ParallelConfTree extends ConfTree {
	
	private class AStarWorker extends Worker {
		
		private int startIndex;
		private int stopIndex;

		public AStarWorker(WorkCrew<AStarWorker> crew) {
			super(crew);
		}

		@Override
		protected void workIt() {
			for (int i=startIndex; i<=stopIndex; i++) {
				scores.set(i, scoreConf(conformations.get(i)));
			}
		}
	}
	
	private WorkCrew<AStarWorker> crew;
	private ArrayList<int[]> conformations;
	private ArrayList<Double> scores;

	public ParallelConfTree(SearchProblem sp, int numThreads) {
		this(sp, sp.pruneMat, sp.useEPIC, numThreads);
	}

	public ParallelConfTree(SearchProblem sp, PruningMatrix pruneMat, boolean useEPIC, int numThreads) {
		super(sp, pruneMat, useEPIC);
		
		// start the work crew
		crew = new WorkCrew<AStarWorker>("AStar");
		for (int i=0; i<numThreads; i++) {
			new AStarWorker(crew);
		}
		crew.start();
		
		conformations = new ArrayList<>();
		scores = new ArrayList<>();
	}
	
	public void askWorkCrewToStop() {
		crew.askToStop();
	}
	
	@Override
	public ArrayList<AStarNode> getChildren(AStarNode curNode) {
		
		if (isFullyAssigned(curNode)) {
			throw new Error("can't expand a fully assigned A* node");
		}
		
		if(curNode.score == Double.POSITIVE_INFINITY) {
			// node impossible, so no children
			return new ArrayList<>();
		}
		
		int nextLevel = nextLevelToExpand(curNode.nodeAssignments);
		
		// clone the conformations and apply the RC changes
		conformations.clear();
		int firstIndex = nextLevel*maxNumRCs;
		for (int i=0; i<maxNumRCs; i++) {
			int rc = unprunedRCsAtPos[firstIndex + i];
			if (rc < 0) {
				break;
			}
			
			int[] childConf = curNode.nodeAssignments.clone();
			childConf[nextLevel] = rc;
			conformations.add(childConf);
		}
		
		// partition conformations among workers
		int numWorkers = crew.getWorkers().size();
		int numConfs = conformations.size();
		int width = (numConfs + numWorkers - 1)/numWorkers;
		int startIndex = 0;
		int stopIndex = startIndex + width - 1;
		for (AStarWorker worker : crew.getWorkers()) {
			worker.startIndex = startIndex;
			worker.stopIndex = stopIndex;
			startIndex += width;
			stopIndex = Math.min(stopIndex + width, numConfs - 1);
		}
		
		// make space for the scores
		scores.clear();
		for (int i=0; i<numConfs; i++) {
			scores.add(null);
		}
		
		// score the conformations
		crew.sendWork();
		try {
			boolean finished = crew.waitForResults(10000);
			if (!finished) {
				throw new Error("Timed out waiting 10 seconds for conformation scoring to finish!"
					+ "\nConformation scoring shouldn't take more than 10 seconds, right?");
			}
		} catch (InterruptedException ex) {
			// something wanted us to stop, so stop, then forward the exception
			throw new Error(ex);
		}
		
		// build the child nodes
		ArrayList<AStarNode> children = new ArrayList<>();
		for (int i=0; i<numConfs; i++) {
			children.add(new AStarNode(conformations.get(i), scores.get(i), useRefinement));
		}
		
		conformations.clear();
		scores.clear();
		
		return children;
	}
}
