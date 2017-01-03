package edu.duke.cs.osprey.minimization;

public interface LineSearcher {
	
	void init(ObjectiveFunction.OneDof f);
	double search(double xd);
	
	// NOTE: when implementing search(), make sure to leave the protein at the pose whose value you return
	
	public static interface NeedsCleanup extends LineSearcher {
		void cleanup();
	}
}
