package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */
public interface KStarScore {

	public enum KStarScoreType {
	    Minimized,//i.e. minimization
	    PairWiseMinimized,//pw min numerator and denominator
	    MinimizedLowerBound,//discrete numerator, pw min denominator
	    MinimizedUpperBound,//pw min numerator, discrete denominator
	    Discrete,//discrete
	    DiscreteLowerBound,//discrete numerator and denominator
	    DiscreteUpperBound;//discrete numerator and denominator
	}
	
	public enum PartitionFunctionType {
		Minimized,//i.e. minimization
		Discrete,//no min; either discrete or pw min
		UpperBound;//1+epsilon on pw min
	}
	
	public MSKStarSettings getSettings();
	public BigDecimal getScore();
	public BigDecimal getLowerBoundScore();
	public BigDecimal getUpperBoundScore();
	
	public String toString();
	public void compute(int maxNumConfs);
	public boolean constrSatisfied();
	public boolean isFullyAssigned();
	public boolean isFinal();
	public boolean isFullyProcessed();
	
}
