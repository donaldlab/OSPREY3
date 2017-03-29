package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;

public interface ResidueOrder {

	public class AAScore {
		public int residuePos;
		public int AATypePos;
		public double score;

		public AAScore(int residuePos, int AATypePos, double score) {
			this.residuePos = residuePos;
			this.AATypePos = AATypePos;
			this.score = score;
		}
	}
	
	public ArrayList<ArrayList<ArrayList<AAScore>>> getNextAssignments(MSSearchProblem[][] objFcnSearch, int numMaxMut);
	
}
