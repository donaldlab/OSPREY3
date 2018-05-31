package edu.duke.cs.osprey.multistatekstar;

import java.io.Serializable;
import java.math.BigDecimal;
import java.util.ArrayList;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public abstract class ResidueOrder implements Serializable {

	public class AAAssignment implements Serializable {
		public int residuePos;
		public int AATypePos;

		public AAAssignment(int residuePos, int AATypePos) {
			this.residuePos = residuePos;
			this.AATypePos = AATypePos;
		}
	}

	public class ResidueAssignment implements Serializable {
		public ArrayList<ArrayList<Integer>> assignments;
		private String string;

		public ResidueAssignment(ArrayList<ArrayList<Integer>> assignments) {
			this.assignments = assignments;
			this.string = null;
		}

		public int length() {
			return assignments.size();
		}

		public ArrayList<Integer> get(int subState) {
			return assignments.get(subState);
		}

		public String toString() {
			if(string!=null) return string;
			
			StringBuilder sb = new StringBuilder();
			for(ArrayList<Integer> al : assignments) {
				for(Integer i : al) sb.append(i + " ");
				sb.append(", ");
			}
			
			string = sb.toString();
			return string;
		}
	}

	public class ResidueAssignmentScore implements Serializable {
		ResidueAssignment assignment;
		BigDecimal score;

		public ResidueAssignmentScore(ResidueAssignment assignment, BigDecimal score) {
			this.assignment = assignment;
			this.score = score;
		}
	}

	public ResidueOrder() {}

	public abstract ArrayList<ArrayList<ArrayList<AAAssignment>>> getNextAssignments(
			LMB objFcn, 
			MSSearchProblem[][] objFcnSearch,
			KStarScore[] objFcnScores,
			int numMaxMut
			);

	public abstract ArrayList<ResidueAssignmentScore> scoreResidueAssignments(
			LMB objFcn, 
			MSSearchProblem[][] objFcnSearch,
			KStarScore[] objFcnScores, 
			int numMaxMut
			);

	protected ArrayList<ArrayList<AAAssignment>> getUnboundAAAssignments(MSSearchProblem bound, 
			ArrayList<ArrayList<AAAssignment>> assignments, MSSearchProblem unbound) {

		ArrayList<ArrayList<AAAssignment>> ans = new ArrayList<>();		
		for(ArrayList<AAAssignment> assignment : assignments) {
			ans.add(new ArrayList<>());
			for(AAAssignment aa : assignment) {//map to unbound state
				int unboundPos=unbound.flexRes.indexOf(bound.flexRes.get(aa.residuePos));
				if(unboundPos != -1)
					ans.get(ans.size()-1).add(new AAAssignment(unboundPos, aa.AATypePos));
			}
			ans.get(ans.size()-1).trimToSize();
		}

		ans.trimToSize();
		return ans;
	}

	protected void getBoundAAAssignmentsHelper(ArrayList<ArrayList<String>> AATypeOptions,
			ArrayList<ArrayList<AAAssignment>> output, ArrayList<Integer> splitPos, 
			String[] wt, String[] buf, int depth, int numMut, int numMaxMut) {

		if(depth==splitPos.size()) {
			ArrayList<AAAssignment> assignment = new ArrayList<>();
			for(int i=0;i<depth;++i) {
				int residuePos = splitPos.get(i);
				int AATypePos = AATypeOptions.get(residuePos).indexOf(buf[i]);
				if(AATypePos == -1)
					throw new RuntimeException("ERROR: AATypeOptions must contain AA");
				assignment.add(new AAAssignment(residuePos, AATypePos));
			}
			assignment.trimToSize();
			output.add(assignment);
			return;
		}

		int residuePos = splitPos.get(depth);
		for(int AATypePos = 0;AATypePos<AATypeOptions.get(residuePos).size();++AATypePos) {
			buf[depth] = AATypeOptions.get(residuePos).get(AATypePos);
			int count = buf[depth].equalsIgnoreCase(wt[residuePos]) ? numMut : numMut+1;
			if(count > numMaxMut) continue;
			getBoundAAAssignmentsHelper(AATypeOptions, output, splitPos, wt, buf, depth+1, count, numMaxMut);
		}
	}
}
