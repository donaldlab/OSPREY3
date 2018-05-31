/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ewakstar;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 *
 * @author Anna Lowegard(anna.lowegard@duke.edu)
 * Adegoke Ojewole (ao68@duke.edu)
 */
public class EWAKScores {

	private ArrayList<EWAKScore> scores;

	public EWAKScores(PartitionFuncDict complexDict, ArrayList<PartitionFuncDict> strandDicts) {

		scores = new ArrayList<>();

		for(String complexSeq : complexDict.getSequences()) {

			EWAKScore score = new EWAKScore();

			for(PartitionFuncDict strandDict : strandDicts) {
				String strandSeq = strandDict.getMatchingSequences(complexSeq).get(0);
				score.add(strandDict.getPartitionFunction(strandSeq));
			}

			score.add(complexDict.getPartitionFunction(complexSeq));
			scores.add(score);
		}

		scores.trimToSize();
	}

	public void sort() {
		Collections.sort(scores, new Comparator<EWAKScore>() {
			@Override
			public int compare(EWAKScore s1, EWAKScore s2) {
				return s1.getScoreLog10().compareTo(s2.getScoreLog10()) <= 0 ? 1 : -1;
			}
		});
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for(EWAKScore score : scores) {
			sb.append(score.toString()+"\n");
		}
		return sb.toString().trim();
	}

}
