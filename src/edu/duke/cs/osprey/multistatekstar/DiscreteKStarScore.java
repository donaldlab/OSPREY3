package edu.duke.cs.osprey.multistatekstar;

public class DiscreteKStarScore implements KStarScore {

	@Override
	public double getScore() {
		// TODO Auto-generated method stub
		return 0;
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
	public double getLowerBoundScore() {
		return getScore();
	}

	@Override
	public double getUpperBoundScore() {
		return getScore();
	}

}
