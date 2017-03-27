package edu.duke.cs.osprey.multistatekstar;

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
	
	public AAScore[][] getNextAssignment(MSSearchProblem[][] objFcnSearch);
	
}
