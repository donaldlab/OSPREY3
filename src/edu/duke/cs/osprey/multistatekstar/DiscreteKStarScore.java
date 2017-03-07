package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;

public class DiscreteKStarScore implements KStarScore {

	@Override
	public BigDecimal getScore() {
		// TODO Auto-generated method stub
		return BigDecimal.ZERO;
	}
	
	public String toString() {
		return null;
	}

	@Override
	public void compute(int maxNumConfs) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public boolean constrSatisfied() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean computed() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public BigDecimal getLowerBoundScore() {
		return getScore();
	}

	@Override
	public BigDecimal getUpperBoundScore() {
		return getScore();
	}

}
