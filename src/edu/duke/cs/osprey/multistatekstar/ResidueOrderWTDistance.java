package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

@SuppressWarnings("serial")
public class ResidueOrderWTDistance extends ResidueOrderStaticMinDomain {

	public static boolean DEBUG = false;
	private HashMap<Integer, HashSet<String>> category2AASet;
	private MSSearchProblem[][] objFcnSearch;
	
	public ResidueOrderWTDistance(MSSearchProblem[][] objFcnSearch) {
		super(objFcnSearch, "product");
		this.objFcnSearch = objFcnSearch;
		this.category2AASet = populateAASet2Category();
	}

	private HashMap<Integer, HashSet<String>> populateAASet2Category() {
		HashMap<Integer, HashSet<String>> ans = new HashMap<>();
		
		HashSet<String> positive = new HashSet<>();
		positive.add("ARG");
		positive.add("HID");
		positive.add("HIE");
		positive.add("HIS");
		positive.add("LYS");
		
		HashSet<String> negative = new HashSet<>();
		negative.add("ASP");
		negative.add("GLU");
		
		HashSet<String> polarUncharged = new HashSet<>();
		polarUncharged.add("SER");
		polarUncharged.add("THR");
		polarUncharged.add("ASN");
		polarUncharged.add("GLN");
		
		HashSet<String> specialCases = new HashSet<>();
		specialCases.add("CYS");
		specialCases.add("SEC");
		specialCases.add("GLY");
		specialCases.add("PRO");
		
		HashSet<String> hydrophobic = new HashSet<>();
		hydrophobic.add("ALA");
		hydrophobic.add("VAL");
		hydrophobic.add("ILE");
		hydrophobic.add("LEU");
		hydrophobic.add("MET");
		hydrophobic.add("PHE");
		hydrophobic.add("TYR");
		hydrophobic.add("TRP");

		int category=0;
		ans.put(category++, positive);
		ans.put(category++, negative);
		ans.put(category++, polarUncharged);
		ans.put(category++, specialCases);
		ans.put(category++, hydrophobic);
		
		return ans;
	}
	
	private int getCategory(String AAType) {
		for(int category : category2AASet.keySet()) {
			if(category2AASet.get(category).contains(AAType))
				return category;
		}
		return -1;
	}
	
	private int getWTDistance(int state, AAAssignment aaAssignment) {
		int complex = objFcnSearch[state].length-1;
		MSSearchProblem search = objFcnSearch[state][complex];
		
		String AAType = search.settings.AATypeOptions.get(aaAssignment.residuePos).get(aaAssignment.AATypePos);
		String WT = search.settings.AATypeOptions.get(aaAssignment.residuePos).get(0);
		
		if(AAType.equals(WT))
			return 0;
		
		int aaCategory = getCategory(AAType);
		int wtCategory = getCategory(WT);
		
		return aaCategory == wtCategory ? 1 : 4;
	}
	
	//return average distance from wt
	protected BigDecimal getBoundStateWTDistance(ResidueAssignment residueAssignment, 
			ArrayList<ArrayList<ArrayList<AAAssignment>>> aaAssignments) {
		//only care about bound state  aa assignments
		int complex = residueAssignment.length()-1;
		int numSplits = aaAssignments.get(0).size();
		double wtDistance = 0;
		
		for(int split=0;split<numSplits;++split) {
			ArrayList<AAAssignment> assignments = aaAssignments.get(complex).get(split);
			for(AAAssignment aaa : assignments) {
				for(int state=0; state<objFcnSearch.length; ++state) {
					wtDistance += getWTDistance(state, aaa);
				}
			}
		}
		
		//avgWTDistance /= (double)numSplits;	
		return BigDecimal.valueOf(wtDistance);
	}
	
	public BigDecimal getResidueAssignmentScore(ResidueAssignment residueAssignment,
			MSSearchProblem[][] objFcnSearch, 
			int numMaxMut) {

		//each residue assignment corresponds to one or more allowed AA assignments.
		//score is based on allowed AA assignments only
		ArrayList<ArrayList<ArrayList<AAAssignment>>> aaAssignments = getAllowedAAAsignments(objFcnSearch, residueAssignment, numMaxMut);

		//store here, so we can retreive best from here later without recomputation
		storeResidue2AAAssignments(residueAssignment, aaAssignments);

		return getBoundStateWTDistance(residueAssignment, aaAssignments);
	}

	@Override
	public ArrayList<ArrayList<ArrayList<AAAssignment>>> getNextAssignments(LMB objFcn,
			MSSearchProblem[][] objFcnSearch, KStarScore[] objFcnScores, int numMaxMut) {

		return super.getNextAssignments(objFcn, objFcnSearch, objFcnScores, numMaxMut);

	}

}
