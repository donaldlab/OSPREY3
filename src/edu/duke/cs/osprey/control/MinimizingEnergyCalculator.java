package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.minimization.ConfMinimizer;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;

public class MinimizingEnergyCalculator implements ConfEnergyCalculator.Async {
	
	// TODO: this should eventually go into a CFP-only area
	// it can be moved when we start refactoring config stuff to prepare for Python-land
	public static MinimizingEnergyCalculator makeFromConfig(SearchProblem search, ConfigFileParser cfp, int queueFactor) {
		int numThreads = cfp.getParams().getInt("MinimizationThreads");
		int numGpus = cfp.getParams().getInt("MinimizationGpus");
		return make(search, numGpus, numThreads, queueFactor);
	}
	
	public static MinimizingEnergyCalculator make(SearchProblem search, int numGpus, int numThreads, int queueFactor) {
		
		int numTasks;
		Factory<? extends EnergyFunction,Molecule> efuncs;
		
		if (numGpus > 0) {
			
			// use gpu-calculated energy functions
			GpuQueuePool gpuPool = new GpuQueuePool(numGpus, 1);
			final GpuEnergyFunctionGenerator egen = new GpuEnergyFunctionGenerator(EnvironmentVars.curEFcnGenerator.ffParams, gpuPool);
			
			efuncs = new Factory<EnergyFunction,Molecule>() {
				@Override
				public EnergyFunction make(Molecule mol) {
					return egen.fullConfEnergy(search.confSpace, search.shellResidues, mol);
				}
			};
			
			numTasks = gpuPool.getNumQueues();
			
		} else {
			
			// plain ol' cpu-calculated energy functions
			efuncs = new Factory<EnergyFunction,Molecule>() {
				@Override
				public EnergyFunction make(Molecule mol) {
					return EnvironmentVars.curEFcnGenerator.fullConfEnergy(search.confSpace, search.shellResidues, mol);
				}
			};
		
			numTasks = numThreads;
		}
		
		// make the thread pool
		TaskExecutor tasks;
		if (numTasks == 0) {
			tasks = new TaskExecutor();
		} else {
			ThreadPoolTaskExecutor poolTasks = new ThreadPoolTaskExecutor();
			poolTasks.start(numTasks, queueFactor);
			tasks = poolTasks;
		}
		
		return new MinimizingEnergyCalculator(search, efuncs, tasks, true);
	}
		
	private SearchProblem search;
	private TaskExecutor tasks;
	private boolean cleanupTasks;
	private ConfMinimizer.Async minimizer;
	
	public MinimizingEnergyCalculator(SearchProblem search, Factory<? extends EnergyFunction,Molecule> efuncs) {
		this(search, efuncs, new TaskExecutor(), true);
	}
	
	public MinimizingEnergyCalculator(SearchProblem search, Factory<? extends EnergyFunction,Molecule> efuncs, TaskExecutor tasks, boolean cleanupTasks) {
		this.search = search;
		this.tasks = tasks;
		this.cleanupTasks = cleanupTasks;
		
		minimizer = new ConfMinimizer.Async(efuncs, search.confSpace, tasks);
	}
	
	@Override
	public int getParallelism() {
		return tasks.getParallelism();
	}
	
	private EnergiedConf postProcessConf(EnergiedConf econf) {
		
		// add post-minimization energy modifications
		if (search.useERef) {
			econf.offsetEnergy(-search.emat.geteRefMat().confERef(econf.getAssignments()));
		}
		if (search.addResEntropy) {
			econf.offsetEnergy(search.confSpace.getConfResEntropy(econf.getAssignments()));
		}
		
		return econf;
	}

	@Override
	public EnergiedConf calcEnergy(ScoredConf conf) {
		return postProcessConf(minimizer.minimizeSync(conf));
	}
	
	@Override
	public void calcEnergyAsync(ScoredConf conf, Listener listener) {
		minimizer.minimizeAsync(conf, new ConfMinimizer.Async.Listener() {
			@Override
			public void onMinimized(EnergiedConf econf) {
				listener.onEnergy(postProcessConf(econf));
			}
		});
	}
	
	@Override
	public void waitForSpace() {
		minimizer.waitForSpace();
	}
	
	@Override
	public void waitForFinish() {
		minimizer.waitForFinish();
	}

	@Override
	public void cleanup() {
		minimizer.cleanup();
		
		if (cleanupTasks && tasks instanceof TaskExecutor.NeedsCleanup) {
			((TaskExecutor.NeedsCleanup)tasks).cleanup();
		}
	}
}
