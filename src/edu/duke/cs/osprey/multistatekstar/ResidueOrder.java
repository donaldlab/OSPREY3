package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.util.ArrayList;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public abstract class ResidueOrder {

	public class AAAssignment {
		public int residuePos;
		public int AATypePos;

		public AAAssignment(int residuePos, int AATypePos) {
			this.residuePos = residuePos;
			this.AATypePos = AATypePos;
		}
	}
	
	public class ResidueAssignment {
		public ArrayList<ArrayList<Integer>> assignments;
		
		public ResidueAssignment(ArrayList<ArrayList<Integer>> assignments) {
			this.assignments = assignments;
		}
		
		public int length() {
			return assignments.size();
		}
		
		public ArrayList<Integer> get(int subState) {
			return assignments.get(subState);
		}
	}
	
	public class ResidueAssignmentScore {
		ResidueAssignment assignment;
		BigDecimal score;
		
		public ResidueAssignmentScore(ResidueAssignment assignment, BigDecimal score) {
			this.assignment = assignment;
			this.score = score;
		}
	}
	
	public class ResisueOrderWorker {
		public int state;
		public int subState;
		public MSSearchProblem search;
		
		public ResisueOrderWorker(MSSearchProblem search, int state, int subState) {
			this.search = search;
			this.state = state;
			this.subState = subState;
		}
	}
	
	public ResidueOrder() {}
	
	public abstract ArrayList<ArrayList<ArrayList<AAAssignment>>> getNextResidueAssignment(
			LMB objFcn, 
			MSSearchProblem[][] objFcnSearch, 
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
