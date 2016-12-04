package edu.duke.cs.osprey.minimization;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.minimization.ConfMinimizer.Async;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.Progress;

public abstract class SpecializedConfMinimizer {
	
	private ThreadPoolTaskExecutor tasks;
	private ConfMinimizer.Async asyncMinimizer;
	private boolean reportProgress;
	
	protected SpecializedConfMinimizer() {
		tasks = null;
		asyncMinimizer = null;
		reportProgress = false;
	}
	
	protected void init(int numThreads, Factory<? extends EnergyFunction,Molecule> efuncs, Factory<? extends Minimizer,MoleculeModifierAndScorer> minimizers, ConfSpace confSpace) {
		
		// start the thread pool
		tasks = new ThreadPoolTaskExecutor();
		tasks.start(numThreads);
		
		// make the minimizer
		asyncMinimizer = new ConfMinimizer.Async(efuncs, confSpace, tasks, minimizers);
	}
	
	public void setReportProgress(boolean val) {
		reportProgress = val;
	}
	
	/**
	 * NOTE: don't call this in a loop, you'll loose all the parallelism
	 * but it's here if you need one-off minimizations
	 */
	public EnergiedConf minimize(ScoredConf conf) {
		return minimize(Arrays.asList(conf)).get(0);
	}
	
	public List<EnergiedConf> minimize(List<ScoredConf> confs) {
		
		// allocate space to hold the minimized values
		List<EnergiedConf> econfs = new ArrayList<>(confs.size());
		for (int i=0; i<confs.size(); i++) {
			econfs.add(null);
		}
		
		// track progress if desired
		final Progress progress;
		if (reportProgress) {
			progress = new Progress(confs.size());
		} else {
			progress = null;
		}
		
		// minimize them all
		for (int i=0; i<confs.size(); i++) {
			final int fi = i;
			asyncMinimizer.minimizeAsync(confs.get(i), new Async.Listener() {
				@Override
				public void onMinimized(EnergiedConf econf) {
					econfs.set(fi, econf);
					if (progress != null) {
						progress.incrementProgress();
					}
				}
			});
		}
		asyncMinimizer.waitForFinish();
		
		return econfs;
	}
	
	public void cleanup() {
		tasks.stop();
		asyncMinimizer.cleanup();
	}
}
