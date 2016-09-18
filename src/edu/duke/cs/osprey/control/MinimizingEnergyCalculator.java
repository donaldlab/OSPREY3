package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.minimization.ConfMinimizer;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;

public class MinimizingEnergyCalculator implements ConfEnergyCalculator.Async {
		
	private SearchProblem search;
	private TaskExecutor tasks;
	private ConfMinimizer.Async minimizer;
	
	public MinimizingEnergyCalculator(SearchProblem search, Factory<? extends EnergyFunction,Molecule> efuncs) {
		this(search, efuncs, 0);
	}
	
	public MinimizingEnergyCalculator(SearchProblem search, Factory<? extends EnergyFunction,Molecule> efuncs, int numThreads) {
		this.search = search;
		
		// make the task executor
		if (numThreads == 0) {
			
			// use the current thread
			tasks = new TaskExecutor();
			
		} else {
			
			// make a thread pool
			ThreadPoolTaskExecutor threadPoolTasks = new ThreadPoolTaskExecutor();
			threadPoolTasks.start(numThreads);
			tasks = threadPoolTasks;
		}
		
		minimizer = new ConfMinimizer.Async(efuncs, search.confSpace, tasks);
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
	public void setListener(Listener listener) {
		minimizer.setListener(new ConfMinimizer.Async.Listener() {
			@Override
			public void onMinimized(EnergiedConf econf, Integer id) {
				listener.onEnergy(postProcessConf(econf));
			}
		});
	}
	
	@Override
	public void calcEnergyAsync(ScoredConf conf) {
		minimizer.minimizeAsync(conf);
	}
	
	@Override
	public void waitForFinish() {
		minimizer.waitForFinish();
	}

	@Override
	public void cleanup() {
		minimizer.cleanup();
		
		if (tasks instanceof ThreadPoolTaskExecutor) {
			((ThreadPoolTaskExecutor)tasks).stop();
		}
	}
}
