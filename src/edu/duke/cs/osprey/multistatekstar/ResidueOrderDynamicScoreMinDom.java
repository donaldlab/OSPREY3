package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

import edu.duke.cs.osprey.tools.ObjectIO;

@SuppressWarnings("serial")
public class ResidueOrderDynamicScoreMinDom extends ResidueOrder {

	public static boolean DEBUG = false;
	
	private ResidueOrderDynamicScore dynamic;
	private ResidueOrderStaticMinDomain mindom;
	//private ResidueOrderWTDistance wtdist;
	
	private double dCoeff = 4.0;
	private double mCoeff = 1.0;
	//private double wCoeff = 1.0;//coeff for dist from wt
	
	private HashMap<ResidueAssignment, Integer> d2A2Rank;
	private HashMap<ResidueAssignment, Integer> m2A2Rank;
	//private HashMap<ResidueAssignment, Integer> w2A2Rank;
	
	public ResidueOrderDynamicScoreMinDom(MSSearchProblem[][] objFcnSearch, 
			ResidueOrderDynamicScore dynamic, double dCoeff, 
			ResidueOrderStaticMinDomain mindom, double mCoeff) {
		super();
		
		this.dCoeff = dCoeff;
		this.mCoeff = mCoeff;
		
		this.dynamic = dynamic;
		this.mindom = mindom;
		//this.wtdist = new ResidueOrderWTDistance(objFcnSearch);
		
		this.d2A2Rank = new HashMap<>();
		this.m2A2Rank = new HashMap<>();
		//this.w2A2Rank = new HashMap<>();
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
		mindom.clearResidue2AAAssignments();
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
		
		/*
		@SuppressWarnings("unchecked")
		ArrayList<ResidueAssignmentScore> w = (ArrayList<ResidueAssignmentScore>) ObjectIO.deepCopy(d);
		wtdist.clearResidue2AAAssignments();
		//set scores to domain sizes
		for(ResidueAssignmentScore ras : w) {
			ras.score = wtdist.getResidueAssignmentScore(ras.assignment, objFcnSearch, numMaxMut);
		}
		
		ArrayList<ResidueAssignmentScore> w2 = new ArrayList<>(w);
		Collections.sort(w2, new Comparator<ResidueAssignmentScore>() {
			@Override
			public int compare(ResidueAssignmentScore ra1, ResidueAssignmentScore ra2) {
				return ra1.score.compareTo(ra2.score)<=0 ? -1 : 1;
			}
		});
		*/
		
		//store ranks; d2, m2, w2 elements dont correspond
		for(int i=0;i<d2.size();++i) {			
			d2A2Rank.put(d2.get(i).assignment, i);
			m2A2Rank.put(m2.get(i).assignment, i);
			//w2A2Rank.put(w2.get(i).assignment, i);
		}
		
		double maxScore = Double.NEGATIVE_INFINITY;
		int maxPos = 0;
		
		for(int i=0;i<d.size();++i) {
			double score = dCoeff * d2A2Rank.get(d.get(i).assignment) 
					+ mCoeff * m2A2Rank.get(m.get(i).assignment);
					//+ wCoeff * w2A2Rank.get(w.get(i).assignment);
			
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
		//w2A2Rank.clear();
		
		return d.get(maxPos).assignment;
	}
	
	@Override
	public ArrayList<ArrayList<ArrayList<AAAssignment>>> getNextAssignments(
			LMB objFcn,
			MSSearchProblem[][] objFcnSearch, 
			KStarScore[] objFcnScores, 
			int numMaxMut
			) {
		ArrayList<ResidueAssignmentScore> dynScores = dynamic.scoreAllResidueAssignments(objFcn, objFcnSearch, objFcnScores, numMaxMut);	
		
		//order unassigned residues by score and return best residue
		ResidueAssignment best = getBestResidueAssignment(dynScores, objFcnSearch, numMaxMut);
		
		//convert to aas that don't violate the allowed number of mutations
		//could also call mindom.getBestAAAssignments
		return dynamic.getAllowedAAAsignments(objFcnSearch, best, numMaxMut);
	}

}
