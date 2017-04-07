package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskListener;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;

public interface FragmentEnergyCalculator {
	
	public static interface InteractionsFactory extends Factory<ForcefieldInteractions,Molecule> {
		// nothing else to do
	}
	
	double calcEnergy(RCTuple frag, InteractionsFactory intersFactory);
	
	// use asynchronous techniques so we can parallelize fragment evaluation
	public static interface Async extends FragmentEnergyCalculator {
		
		public static interface Listener extends TaskListener<Double> {
			// nothing else to do
		}
		
		void calcEnergyAsync(RCTuple frag, InteractionsFactory intersFactory, Listener listener);
		TaskExecutor getTasks();
		void cleanup();
	}
}
