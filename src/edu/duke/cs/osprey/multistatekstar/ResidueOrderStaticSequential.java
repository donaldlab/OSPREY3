/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */

public class ResidueOrderStaticSequential implements ResidueOrder {

	protected ArrayList<ArrayList<ArrayList<ArrayList<AAScore>>>> scores;

	public ResidueOrderStaticSequential(MSSearchProblem[][] objFcnSearch) {
		scores = null;
		allocate(objFcnSearch);
		init(objFcnSearch);
	}

	private void allocate(MSSearchProblem[][] objFcnSearch) {
		//create a priority queue for each state and substate
		scores = new ArrayList<>();
		for(int state=0;state<objFcnSearch.length;++state) {
			scores.add(new ArrayList<>());
			for(int subState=0;subState<objFcnSearch[state].length;++subState) {
				scores.get(state).add(new ArrayList<>());
				for(int residuePos=0;residuePos<objFcnSearch[state][subState].settings.AATypeOptions.size();++residuePos) {
					scores.get(state).get(subState).add(new ArrayList<>());
				}
				scores.get(state).get(subState).trimToSize();
			}
			scores.get(state).trimToSize();
		}
		scores.trimToSize();
	}

	/**
	 * iterate through all aatypeoptions and assign an increasing numerical score
	 */
	protected void init(MSSearchProblem[][] objFcnSearch) {
		for(int state=0;state<scores.size();++state) {
			int score = -1;
			for(int subState=0;subState<scores.get(state).size();++subState) {
				ArrayList<ArrayList<String>> AATypeOptions = objFcnSearch[state][subState].settings.AATypeOptions;
				for(int residuePos=0;residuePos<AATypeOptions.size();++residuePos) {
					ArrayList<String> AATypes = AATypeOptions.get(residuePos);
					for(int AATypePos=0;AATypePos<AATypes.size();++AATypePos) {
						scores.get(state).get(subState).get(residuePos).add(new AAScore(residuePos, AATypePos, ++score));
					}
				}
			}
		}
	}

	@Override
	public ArrayList<ArrayList<ArrayList<AAScore>>> getNextAssignments(MSSearchProblem[][] objFcnSearch, int numMaxMut) {
		ArrayList<ArrayList<ArrayList<AAScore>>> ans = new ArrayList<>();
		
		int state = 0;
		int numSubStates = objFcnSearch[state].length;
		//bound state is the sequence
		MSSearchProblem bound = objFcnSearch[state][numSubStates-1];

		if(bound.getNumAssignedPos()==0) {//root node; add all allowed single mutations from unbound states
			ArrayList<ArrayList<AAScore>> boundAssignments = getBoundStateAssignments(state, objFcnSearch[state], 0, numMaxMut);
			for(int subState=0;subState<numSubStates-1;++subState) {
				MSSearchProblem unbound = objFcnSearch[state][subState];
				ans.add(getUnboundStateAssignments(bound, boundAssignments, unbound));
			}
			ans.add(boundAssignments);
		}

		else {//add all allowed mutations at the next numerical splitPos
			ArrayList<Integer> splitPos = bound.getPosNums(false);
			if(splitPos.size()==0)
				throw new RuntimeException("ERROR: there are no positions to split");
		}
		
		ans.trimToSize();
		return ans;
	}

	protected ArrayList<ArrayList<AAScore>> nextSplitPosAssignments(int splitPos) {
		return null;
	}

	private ArrayList<ArrayList<AAScore>> getUnboundStateAssignments(MSSearchProblem bound, 
			ArrayList<ArrayList<AAScore>> assignments, MSSearchProblem unbound) {

		ArrayList<ArrayList<AAScore>> ans = new ArrayList<>();		
		for(ArrayList<AAScore> assignment : assignments) {
			ans.add(new ArrayList<>());
			for(AAScore aa : assignment) {//map to unbound state
				int unboundPos=unbound.flexRes.indexOf(bound.flexRes.get(aa.residuePos));
				if(unboundPos != -1)
					ans.get(ans.size()-1).add(new AAScore(unboundPos, aa.AATypePos,-1));
			}
			ans.get(ans.size()-1).trimToSize();
		}

		ans.trimToSize();
		return ans;
	}

	private ArrayList<ArrayList<AAScore>> getBoundStateAssignments(int state, MSSearchProblem[] search, int splitPos, int numMaxMut) {
		ArrayList<Integer> complexPos = new ArrayList<>();
		MSSearchProblem complex = search[search.length-1];
		for(int subState=0;subState<search.length-1;++subState)
			complexPos.add(complex.flexRes.indexOf(search[subState].flexRes.get(splitPos)));
		complexPos.trimToSize();

		ArrayList<ArrayList<AAScore>> ans = new ArrayList<>();
		String[] wt = MSKStarNode.WT_SEQS.get(state);
		String[] buf = new String[wt.length];
		getBoundStateAssignmentsHelper(complex.allowedAAs, ans, complexPos, wt, buf, 0, 0, numMaxMut);

		ans.trimToSize();
		return ans;
	}

	private void getBoundStateAssignmentsHelper(ArrayList<ArrayList<String>> AATypeOptions,
			ArrayList<ArrayList<AAScore>> output, ArrayList<Integer> splitPos, 
			String[] wt, String[] buf, int depth, int numMut, int numMaxMut) {

		if(depth==splitPos.size()) {
			ArrayList<AAScore> assignment = new ArrayList<>();
			for(int i=0;i<depth;++i) {
				int residuePos = splitPos.get(i);
				int AATypePos = AATypeOptions.get(residuePos).indexOf(buf[i]);
				if(AATypePos == -1)
					throw new RuntimeException("ERROR: AATypeOptions must contain AA");
				assignment.add(new AAScore(residuePos, AATypePos, -1));
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
			getBoundStateAssignmentsHelper(AATypeOptions, output, splitPos, wt, buf, depth+1, count, numMaxMut);
		}
	}

}
