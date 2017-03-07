package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;

public interface KStarScore {

	public BigDecimal getScore();
	public BigDecimal getLowerBoundScore();
	public BigDecimal getUpperBoundScore();
	
	public String toString();
	public void compute(int maxNumConfs);
	public boolean constrSatisfied();
	public boolean computed();
	
}
