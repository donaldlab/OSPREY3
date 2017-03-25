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
	private KStarScore[] kss;//kstar score objects for objective function or constraint
	private BigDecimal[] ksvals;//kstar score values
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
		this.ksvals = new BigDecimal[ksvals.length];
		this.tree = tree;
		score = NEGATIVE_INFINITY;
	}

	public BigDecimal getScore() {
		return score;
	}

	public void setScore(BigDecimal val) {
		score = val;
	}

	public BigDecimal[] getKStarValues(LMB constr) {
		kss = getConstrKStarScores(constr);
		for(int i=0;i<kss.length;++i) {
			if(kss[i]==null) ksvals[i] = BigDecimal.ZERO;
			//?
			else ksvals[i] = kss[i].getScore();
		}
		return ksvals;
	}
	
	public KStarScore[] getConstrKStarScores(LMB lmb) {
		BigDecimal[] coeffs = lmb.getCoeffs();

		if(lmb.equals(tree.objFcn)) {//objective function cannot have 0 coefficients
			for(int i=0;i<coeffs.length;++i) {
				if(coeffs[i].compareTo(BigDecimal.ZERO) < 0) kss[i] = kssLB[i];
				else if(coeffs[i].compareTo(BigDecimal.ZERO) > 0) kss[i] = kssUB[i];
				else throw new RuntimeException("ERROR: objective function coefficient cannot be 0");
				if(kss[i] == null) throw new RuntimeException("ERROR: must have either upper or lower bound K* score");
			}
		}

		else {//leave out score when coefficient is 0
			for(int i=0;i<coeffs.length;++i) {
				if(coeffs[i].compareTo(BigDecimal.ZERO) < 0) kss[i] = kssUB[i];
				else if(coeffs[i].compareTo(BigDecimal.ZERO) > 0) kss[i] = kssLB[i];
				else kss[i] = null;
			}
		}
		
		return kss;
	}

	public boolean isLeafNode() {
		kss = getConstrKStarScores(tree.objFcn);
		for(int i=0;i<kss.length;++i) 
			if(!kss[i].isFullyProcessed()) return false;
		return true;
	}
}

