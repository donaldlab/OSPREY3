package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;
import java.util.HashMap;

@SuppressWarnings("serial")
public class ResidueOrderMeta extends ResidueOrder {

	public static boolean DEBUG = false;
	
	private ResidueOrderGMEC gmec;
	private ResidueOrderMinDomain mindom;
	//private ResidueOrderWTDistance wtdist;
	
	private double gCoeff = 1.0;
	private double mCoeff = 1.0;
	
	private HashMap<String, Integer> g2A2Rank;
	private HashMap<String, Integer> m2A2Rank;
	
	public ResidueOrderMeta(ResidueOrderGMEC gmec, double gCoeff, 
			ResidueOrderMinDomain mindom, double mCoeff) {
		super();
		
		//this is intentional! when you say dynamicmindom 4 1, we want the decision
		//to weigh more with the dcoefficient, so we are reducing the magnitude
		//of the mcoefficient.
		this.gCoeff = 1.0/mCoeff;
		this.mCoeff = 1.0/gCoeff;
		
		this.gmec = gmec;
		this.mindom = mindom;
		
		this.g2A2Rank = new HashMap<>();
		this.m2A2Rank = new HashMap<>();
	}
	
	private ResidueAssignment getBestResidueAssignment(
			ArrayList<ResidueAssignmentScore> gmecScores,
			ArrayList<ResidueAssignmentScore> mindomScores
			) {
		
		int size = gmecScores.size();
		if(size != mindomScores.size())
			throw new RuntimeException("ERROR: gmec residue assignment size != mindomain residue assignment size");
		
		//store ranks; d2, m2, w2 elements dont correspond
		for(int i=0;i<size;++i) {			
			g2A2Rank.put(gmecScores.get(i).assignment.toString(), i);
			m2A2Rank.put(mindomScores.get(i).assignment.toString(), i);
		}
		
		double bestMetaRank = Double.POSITIVE_INFINITY;
		int bestPos = 0;
		
		//create new meta rank. lowest rank wins
		for(int i=0;i<size;++i) {
			double metaRank = gCoeff * g2A2Rank.get(gmecScores.get(i).assignment.toString()) 
					+ mCoeff * m2A2Rank.get(gmecScores.get(i).assignment.toString());
			
			if(DEBUG) {
				System.out.println("metaRank: "+metaRank+
						" dCoeff: "+gCoeff+" rank: "+g2A2Rank.get(gmecScores.get(i).assignment.toString())+
						" mCoeff: "+mCoeff+" rank: "+m2A2Rank.get(gmecScores.get(i).assignment.toString()));
			}
			
			if(metaRank<bestMetaRank) {
				bestMetaRank = metaRank;
				bestPos = i;
			}
		}
		
		g2A2Rank.clear();
		m2A2Rank.clear();
		
		return gmecScores.get(bestPos).assignment;
	}
	
	@Override
	public ArrayList<ArrayList<ArrayList<AAAssignment>>> getNextAssignments(
			LMB objFcn,
			MSSearchProblem[][] objFcnSearch, 
			KStarScore[] objFcnScores, 
			int numMaxMut
			) {
		ArrayList<ResidueAssignmentScore> gmecScores = gmec.scoreResidueAssignments(objFcn, objFcnSearch, objFcnScores, numMaxMut);
		ArrayList<ResidueAssignmentScore> mindomScores = mindom.scoreResidueAssignments(objFcn, objFcnSearch, objFcnScores, numMaxMut);
		
		//order unassigned residues by score and return best residue
		ResidueAssignment best = getBestResidueAssignment(gmecScores, mindomScores);
		
		//could equally call mindom.getAAAssignments
		return gmec.getAAAssignments(best);
	}

	@Override
	public ArrayList<ResidueAssignmentScore> scoreResidueAssignments(LMB objFcn, MSSearchProblem[][] objFcnSearch,
			KStarScore[] objFcnScores, int numMaxMut) {
		// TODO Auto-generated method stub
		return null;
	}

}
