package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.util.Arrays;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Node for multistate k* tree
 *
 */
public class MSKStarNode {

	public static final BigDecimal NEGATIVE_INFINITY = new BigDecimal("-1e1024");
	
	private KStarScore[] kssLB;//leaf nodes are stored here
	private KStarScore[] kssUB;//leaf node: array of nulls
	private KStarScore[] kss;//kstar scores for objective function or constraint
	private final MSKStarTree tree;//has all required state
	private BigDecimal score;
	
	public MSKStarNode(
			MSKStarTree tree,
			KStarScore[] kssLB, 
			KStarScore[] kssUB
			) {
		this.kssUB = kssLB;
		this.kssLB = kssUB;
		this.kss = new KStarScore[kssLB.length];
		this.tree = tree;
		score = NEGATIVE_INFINITY;
	}
	
	public BigDecimal getScore() {
		return score;
	}
	
	public void setScore(BigDecimal val) {
		score = val;
	}
	
	public KStarScore[] getObjFcnKStarScores() {
		BigDecimal[] coeffs = tree.objFcn.getCoeffs();
		for(int i=0;i<coeffs.length;++i) {
			if(coeffs[i].compareTo(BigDecimal.ZERO) < 0) kss[i] = kssLB[i];
			else if(coeffs[i].compareTo(BigDecimal.ZERO) > 0) kss[i] = kssUB[i];
			else throw new RuntimeException("ERROR: objective function coefficient cannot be 0");
			if(kss[i] == null) throw new RuntimeException("ERROR: must have either upper or lower bound K* score");
		}
		return kss;
	}
	
	public KStarScore[] getConstrKStarScores(LMB constr) {
		BigDecimal[] coeffs = constr.getCoeffs();
		Arrays.fill(kss, null);
		//leave out scores when coefficient is 0
		for(int i=0;i<coeffs.length;++i) {
			if(coeffs[i].compareTo(BigDecimal.ZERO) < 0) kss[i] = kssUB[i];
			else if(coeffs[i].compareTo(BigDecimal.ZERO) > 0) kss[i] = kssLB[i];
		}
		return kss;
	}
	
	public boolean isLeafNode() {
		kss = getObjFcnKStarScores();
		for(int i=0;i<kss.length;++i) 
			if(!kss[i].isFullyProcessed()) return false;
		return true;
	}
}

