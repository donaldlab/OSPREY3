package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskListener;
import edu.duke.cs.osprey.tools.AutoCleanable;

public interface FragmentEnergyCalculator {
	
	SimpleConfSpace getConfSpace();
	double calcEnergy(RCTuple frag, ResidueInteractions inters);
	
	// use asynchronous techniques so we can parallelize fragment evaluation
	public static interface Async extends FragmentEnergyCalculator, AutoCleanable {
		
		public static interface Listener extends TaskListener<Double> {
			// nothing else to do
		}
		
		void calcEnergyAsync(RCTuple frag, ResidueInteractions inters, Listener listener);
		TaskExecutor getTasks();
	}
}
