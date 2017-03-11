package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;

public interface KStarScore {

	public enum KStarScoreType {
	    Continuous,//i.e. minimization
	    Discrete,//discrete
	    DiscretePairWise,//pw min numerator and denominator
	    DiscreteLowerBound,//discrete numerator, pw min denominator
	    DiscreteUpperBound;//pw min numerator, discrete denominator
	}
	
	public enum PartitionFunctionType {
		Continuous,//i.e. minimization
		Discrete,//no min; either discrete or pw min
		DiscreteUpperBound;//1+epsilon on pw min
	}
	
	public BigDecimal getScore();
	public BigDecimal getLowerBoundScore();
	public BigDecimal getUpperBoundScore();
	
	public String toString();
	public void compute(int maxNumConfs);
	public boolean constrSatisfied();
	public boolean computed();
	
}
