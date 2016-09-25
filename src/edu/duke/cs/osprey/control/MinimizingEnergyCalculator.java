package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.minimization.ConfMinimizer;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;

public class MinimizingEnergyCalculator implements ConfEnergyCalculator.Async {
		
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
