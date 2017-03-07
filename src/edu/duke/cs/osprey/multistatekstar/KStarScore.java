package edu.duke.cs.osprey.multistatekstar;

public interface KStarScore {

	public double getScore();
	public double getLowerBoundScore();
	public double getUpperBoundScore();
	
	public String toString();
	public void compute(int maxNumConfs);
	public boolean constrSatisfied();
	public boolean computed();
	
}
