package edu.duke.cs.osprey.ematrix;

import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.structure.Residue;

public class ExampleParallelEmat {
	
	public static EnergyMatrix calcEmat(ConfSpace confSpace, List<Residue> shellResidues, int numThreads) {
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		
		// build the energy matrix calculator
		// (you don't have to worry about concurrency issues here, SimpleEnergyCalculator takes care of all of that)
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, confSpace, shellResidues);
		SimpleEnergyMatrixCalculator ematcalc = new SimpleEnergyMatrixCalculator(ecalc);
		
		// start a thread pool of whatever size you want
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(numThreads);
	
		// calculate the emat
		EnergyMatrix emat = ematcalc.calcEnergyMatrix(tasks);
		
		// cleanup
		tasks.stop();
		
		return emat;
	}
}
