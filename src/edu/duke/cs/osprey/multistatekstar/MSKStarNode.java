package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Node for multistate k* tree
 *
 */
public class MSKStarNode {

	public static final BigDecimal NEGATIVE_INFINITY = new BigDecimal("-1e1024");
	
	private KStarScore[] kssLB;
	private KStarScore[] kssUB;
	private final MSKStarTree tree;//has all required state
	private BigDecimal score;
	
	public MSKStarNode(
			MSKStarTree tree,
			KStarScore[] kssLB, 
			KStarScore[] kssUB
			) {
		this.kssUB = kssLB;
		this.kssLB = kssUB;
		this.tree = tree;
		score = NEGATIVE_INFINITY;
	}
	
	public BigDecimal getScore() {
		return score;
	}
	
	public void setScore(BigDecimal val) {
		score = val;
	}
}

