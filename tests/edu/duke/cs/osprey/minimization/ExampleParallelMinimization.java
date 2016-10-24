package edu.duke.cs.osprey.minimization;

import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Factory;

public class ExampleParallelMinimization {
	
	public static List<EnergiedConf> minimize(ConfSpace confSpace, List<Residue> shellResidues, List<ScoredConf> confs) {
		
		// make the energy function factory
		// (each thread gets its own energy function and copy of the molecule, to prevent racing)
		GpuEnergyFunctionGenerator egen = new GpuEnergyFunctionGenerator(EnvironmentVars.curEFcnGenerator.ffParams, new GpuQueuePool());
		Factory<EnergyFunction,Molecule> efuncs = new Factory<EnergyFunction,Molecule>() {
			@Override
			public EnergyFunction make(Molecule mol) {
				return egen.fullConfEnergy(confSpace, shellResidues, mol);
			}
		};
		
		// make a thread pool to match the gpu pool
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(egen.getOpenclQueuePool().getNumQueues());
		
		List<EnergiedConf> minimizedConfs = new ConfMinimizer().minimize(confs, efuncs, confSpace, tasks);
		
		// cleanup
		egen.cleanup();
		tasks.stop();
		
		return minimizedConfs;
	}
}
