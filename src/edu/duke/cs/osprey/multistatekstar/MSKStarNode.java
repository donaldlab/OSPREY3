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

import java.math.BigDecimal;
import java.util.ArrayList;

import edu.duke.cs.osprey.control.ParamSet;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Node for multistate k* tree
 *
 */
public class MSKStarNode {

	public static LMB OBJ_FUNC;//global objective function that we are optimizing for
	public static ArrayList<String[]> WT_SEQS;//wild-type bound state sequences for all states
	public static int NUM_MAX_MUT;//maximum number of allowed mutations

	private KStarScore[] ksLB;//lower bound k* objects
	private KStarScore[] ksUB;//upper bound k* objects
	private KStarScore[] ksObjFunc;//k* objects that minimize objective function
	private BigDecimal[] kss;//kstar score values
	private BigDecimal score;//objective function value; smaller is better

	public MSKStarNode(
			KStarScore[] ksLB, 
			KStarScore[] ksUB
			) {
		this.ksLB = ksLB;
		this.ksUB = ksUB;
		this.ksObjFunc = new KStarScore[ksLB.length];
		this.kss = new BigDecimal[ksLB.length];
		score = null;
	}
	
	public boolean isFullyAssigned() {
		return ksLB[0].isFullyAssigned();
	}
	
	public ArrayList<MSKStarNode> split(ParamSet msParams) {
		ksObjFunc = getStateKStarObjects(OBJ_FUNC);
		MSSearchProblem[][] objFcnSearch = new MSSearchProblem[ksObjFunc.length][];
		for(int state=0;state<objFcnSearch.length;++state) objFcnSearch[state] = ksObjFunc[state].getSettings().search;
			
		ResidueOrder order = ResidueOrderFactory.getResidueOrder(msParams, objFcnSearch);
		order.getNextAssignments(objFcnSearch, NUM_MAX_MUT-getNumMutations(0));
		return null;
	}
	
	private int getNumMutations(int state) {
		KStarScore score = ksLB[state];
		int boundState = score.getSettings().search.length-1;
		ArrayList<ArrayList<String>> AATypeOptions = score.getSettings().search[boundState].settings.AATypeOptions;
		
		int ans = 0;
		for(int pos : score.getSettings().search[boundState].getPosNums(true)) {
			if(AATypeOptions.get(pos).size()>0) ans++;
			else if(!AATypeOptions.get(pos).get(0).equalsIgnoreCase(WT_SEQS.get(state)[pos]))
				ans++;
		}
		return ans;
	}

	public BigDecimal getScore() {
		return score;
	}
	
	public void setScore(LMB lmb) {
		score = lmb.eval(getStateKStarScores(lmb));
	}

	public BigDecimal[] getStateKStarScores(LMB lmb) {
		BigDecimal[] coeffs = lmb.getCoeffs();
		//always want a lower bound on the lmb value
		for(int i=0;i<coeffs.length;++i) {
			if(coeffs[i].compareTo(BigDecimal.ZERO) < 0) kss[i] = ksUB[i].getUpperBoundScore();
			else if(coeffs[i].compareTo(BigDecimal.ZERO) > 0) kss[i] = ksLB[i].getLowerBoundScore();
			else {//coeffs[i]==0
				if(lmb.equals(OBJ_FUNC)) throw new RuntimeException("ERROR: objective function coefficient cannot be 0");
				else kss[i] = BigDecimal.ZERO;
			}
		}
		return kss;
	}
	
	public KStarScore[] getStateKStarObjects(LMB lmb) {
		BigDecimal[] coeffs = lmb.getCoeffs();
		//always want a lower bound on the lmb value
		for(int i=0;i<coeffs.length;++i) {
			if(coeffs[i].compareTo(BigDecimal.ZERO) < 0) ksObjFunc[i] = ksUB[i];
			else if(coeffs[i].compareTo(BigDecimal.ZERO) > 0) ksObjFunc[i] = ksLB[i];
			else {//coeffs[i]==0
				if(lmb.equals(OBJ_FUNC)) throw new RuntimeException("ERROR: objective function coefficient cannot be 0");
				else ksObjFunc = null;
			}
		}
		return ksObjFunc;
	}

	public boolean isLeafNode() {
		//kssLB is never null for fully defined sequences, minimized or not
		for(int i=0;i<ksLB.length;++i) if(!ksLB[i].isFullyProcessed()) return false;
		return true;
	}
	
	public boolean constrSatisfied() {
		for(KStarScore score : getStateKStarObjects(OBJ_FUNC)) {
			if(score==null) continue;
			if(!score.constrSatisfied()) return false;
		}
		return true;
	}
	
	public String toString() {
		return "";
	}
}
