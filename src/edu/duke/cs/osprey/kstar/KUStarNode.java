package edu.duke.cs.osprey.kstar;

import java.math.BigDecimal;


/*
 * TODO
 * 1) set score
 * 2) when is this node fully assigned?
 * 3) when does this score need refinement?
 */

public class KUStarNode {

	private KSCalc calc = null;
	private BigDecimal score = null;
	
	public KUStarNode() {

	}
	
	BigDecimal getScore() {
		if( score == null ) score = BigDecimal.ONE.divide(calc.getKStarScore(), 4);
		return score;
	}
	
	void expand() {
		
	}
	
	boolean scoreNeedsRefinement() {
		return false;
	}
}
