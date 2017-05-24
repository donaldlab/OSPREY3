package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

import edu.duke.cs.osprey.tools.ObjectIO;

@SuppressWarnings("serial")
public class ResidueOrderDynamicScoreMinDom extends ResidueOrder {

	private ResidueOrderDynamicScore dynamic;
	private ResidueOrderStaticMinDomain mindom;
	
	private double dCoeff = 1.0;
	private double mCoeff = 2.0;
	
	private HashMap<ResidueAssignment, Integer> d2A2Rank;
	private HashMap<ResidueAssignment, Integer> m2A2Rank;
	public static boolean DEBUG = false;
	
	public ResidueOrderDynamicScoreMinDom(MSSearchProblem[][] objFcnSearch, double dCoeff, double mCoeff) {
		super();
		
		this.dCoeff = dCoeff;
		this.mCoeff = mCoeff;
		
		dynamic = new ResidueOrderDynamicScore(objFcnSearch, "discrepancy");
		mindom = new ResidueOrderStaticMinDomain(objFcnSearch, "product");
		
		d2A2Rank = new HashMap<>();
		m2A2Rank = new HashMap<>();
	}
	
	private ResidueAssignment getBestResidueAssignment(
			ArrayList<ResidueAssignmentScore> d,
			MSSearchProblem[][] objFcnSearch,
			int numMaxMut) {	
		ArrayList<ResidueAssignmentScore> d2 = new ArrayList<>(d);
		Collections.sort(d2, new Comparator<ResidueAssignmentScore>() {
			@Override
			public int compare(ResidueAssignmentScore ra1, ResidueAssignmentScore ra2) {
				return ra1.score.compareTo(ra2.score)>=0 ? -1 : 1;
			}
		});
		
		@SuppressWarnings("unchecked")
		ArrayList<ResidueAssignmentScore> m = (ArrayList<ResidueAssignmentScore>) ObjectIO.deepCopy(d);
		//set scores to domain sizes
		for(ResidueAssignmentScore ras : m) {
			ras.score = mindom.getResidueAssignmentScore(ras.assignment, objFcnSearch, numMaxMut);
		}
		
		ArrayList<ResidueAssignmentScore> m2 = new ArrayList<>(m);
		Collections.sort(m2, new Comparator<ResidueAssignmentScore>() {
			@Override
			public int compare(ResidueAssignmentScore ra1, ResidueAssignmentScore ra2) {
				return ra1.score.compareTo(ra2.score)<=0 ? -1 : 1;
			}
		});
		
		for(int i=0;i<d2.size();++i) {			
			d2A2Rank.put(d2.get(i).assignment, i);
			m2A2Rank.put(m2.get(i).assignment, i);
		}
		
		double maxScore = Double.NEGATIVE_INFINITY;
		int maxPos = 0;
		
		for(int i=0;i<d.size();++i) {
			double score = dCoeff * d2A2Rank.get(d.get(i).assignment) + mCoeff * m2A2Rank.get(m.get(i).assignment);
			
			if(DEBUG) {
				System.out.println("score: "+score+" dCoeff: "+dCoeff+" rank: "+d2A2Rank.get(d.get(i).assignment)+" mCoeff: "+mCoeff+" rank: "+m2A2Rank.get(m.get(i).assignment));
			}
			
			if(score>maxScore) {
				maxScore = score;
				maxPos = i;
			}
		}
		
		d2A2Rank.clear();
		m2A2Rank.clear();
		
		return d.get(maxPos).assignment;
	}
	
	@Override
	public ArrayList<ArrayList<ArrayList<AAAssignment>>> getNextAssignments(
			LMB objFcn,
			MSSearchProblem[][] objFcnSearch, 
			KStarScore[] objFcnScores, 
			int numMaxMut
			) {
		ArrayList<ResidueAssignmentScore> dynScores = dynamic.getAllResidueAssignments(objFcn, objFcnSearch, objFcnScores, numMaxMut);	
		
		//order unassigned residues by score and return best residue
		ResidueAssignment best = getBestResidueAssignment(dynScores, objFcnSearch, numMaxMut);
		
		//convert to aas that don't violate the allowed number of mutations
		//could also call mindom.getBestAAAssignments
		return dynamic.getAllowedAAAsignments(objFcnSearch, best, numMaxMut);
	}

}
