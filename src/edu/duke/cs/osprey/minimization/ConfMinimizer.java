package edu.duke.cs.osprey.minimization;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.tools.Profiler;
import edu.duke.cs.osprey.tools.Progress;

public class ConfMinimizer {
	
	private Factory<Minimizer,MoleculeModifierAndScorer> minimizers;
	
	public ConfMinimizer() {
		this(new Factory<Minimizer,MoleculeModifierAndScorer>() {
			@Override
			public Minimizer make(MoleculeModifierAndScorer mof) {
				return new CCDMinimizer(mof, true);
			}
		});
	}
	
	public ConfMinimizer(Factory<Minimizer,MoleculeModifierAndScorer> minimizers) {
		this.minimizers = minimizers;
	}
	
	public Minimizer.Result minimize(ParameterizedMoleculeCopy pmol, int[] conf, EnergyFunction efunc, ConfSpace confSpace) {
		
		// TEMP
		Profiler p = new Profiler();
		p.start("mof");
		
		RCTuple tuple = new RCTuple(conf);
		MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, confSpace, tuple, pmol);
		
		// TEMP
		p.start("minimizer");
		
		Minimizer minimizer = minimizers.make(mof);
		
		// TEMP
		p.start("minimize");
		
		Minimizer.Result result = minimizer.minimize();
		
		// TEMP
		p.start("cleanup minimizer");
		
		if (minimizer instanceof Minimizer.NeedsCleanup) {
			((Minimizer.NeedsCleanup)minimizer).cleanup();
		}
		
		// TEMP
		p.start("cleanup mof");
		
		mof.cleanup();
		
		// TEMP
		p.stop();
		System.out.println("ConfMinimizer " + p.makeReport(TimeUnit.MILLISECONDS));
		
		return result;
	}
	
	public EnergiedConf minimize(ParameterizedMoleculeCopy pmol, ScoredConf conf, EnergyFunction efunc, ConfSpace confSpace) {
		Minimizer.Result result = minimize(pmol, conf.getAssignments(), efunc, confSpace);
		return new EnergiedConf(conf, result.energy);
	}
	
	public List<EnergiedConf> minimize(ParameterizedMoleculeCopy pmol, List<ScoredConf> confs, EnergyFunction efunc, ConfSpace confSpace) {
		Progress progress = new Progress(confs.size());
		List<EnergiedConf> econfs = new ArrayList<>();
		for (ScoredConf conf : confs) {
			econfs.add(minimize(pmol, conf, efunc, confSpace));
			progress.incrementProgress();
		}
		return econfs;
	}
	
	public List<EnergiedConf> minimize(List<ScoredConf> confs, Factory<? extends EnergyFunction,Molecule> efuncs, ConfSpace confSpace) {
		return minimize(confs, efuncs, confSpace, new TaskExecutor());
	}
	
	public List<EnergiedConf> minimize(List<ScoredConf> confs, Factory<? extends EnergyFunction,Molecule> efuncs, ConfSpace confSpace, TaskExecutor tasks) {
		
		final Progress progress = new Progress(confs.size());
		
		List<EnergiedConf> econfs = new ArrayList<>(confs.size());
		for (int i=0; i<confs.size(); i++) {
			econfs.add(null);
		}
		
		Async async = new Async(efuncs, confSpace, tasks, this);
		
		// minimize them all
		for (int i=0; i<confs.size(); i++) {
			final int fi = i;
			async.minimizeAsync(confs.get(i), new Async.Listener() {
				@Override
				public void onMinimized(EnergiedConf econf) {
					econfs.set(fi, econf);
					progress.incrementProgress();
				}
			});
		}
		async.waitForFinish();
		async.cleanup();
		
		return econfs;
	}
	
	public static class Async {
		
		private static class EfuncAndPMC {
			public ParameterizedMoleculeCopy pmol;
			public EnergyFunction efunc;
		}
		
		private class Task implements Runnable {
			
			public ScoredConf conf;
			public EnergiedConf minimizedConf;
			
			@Override
			public void run() {
				minimizedConf = minimizeSync(conf);
			}
		}
		
		public static interface Listener {
			void onMinimized(EnergiedConf econf);
		}
	
		private ConfSpace confSpace;
		private TaskExecutor tasks;
		private ConfMinimizer confMinimizer;
		private ObjectPool<EfuncAndPMC> pool;
	
		public Async(Factory<? extends EnergyFunction,Molecule> efuncs, ConfSpace confSpace, TaskExecutor tasks) {
			this(efuncs, confSpace, tasks, new ConfMinimizer());
		}
		
		public Async(Factory<? extends EnergyFunction,Molecule> efuncs, ConfSpace confSpace, TaskExecutor tasks, ConfMinimizer confMinimizer) {
			
			this.confSpace = confSpace;
			this.tasks = tasks;
			this.confMinimizer = confMinimizer;
			
			// make a pool for molecules and energy functions
			// to keep concurrent tasks from racing each other
			pool = new ObjectPool<>(new Factory<EfuncAndPMC,Void>() {
				@Override
				public EfuncAndPMC make(Void context) {
					
					EfuncAndPMC out = new EfuncAndPMC();
					out.pmol = new ParameterizedMoleculeCopy(confSpace);
					out.efunc = efuncs.make(out.pmol.getCopiedMolecule());
					return out;
				}
			});
		}
		
		public EnergiedConf minimizeSync(ScoredConf conf) {
			EfuncAndPMC em;
			synchronized (pool) {
				em = pool.checkout();
			}
			try {
				return confMinimizer.minimize(em.pmol, conf, em.efunc, confSpace);
			} finally {
				synchronized (pool) {
					pool.release(em);
				}
			}
		}
		
		public void minimizeAsync(ScoredConf conf, Listener listener) {
			
			if (listener == null) {
				throw new IllegalArgumentException("listener can't be null");
			}
			
			// make the task to handle the minimization
			Task task = new Task();
			task.conf = conf;
			
			// submit the task with a chaining task listener
			tasks.submit(task, new TaskExecutor.TaskListener() {
				@Override
				public void onFinished(Runnable taskBase) {
					listener.onMinimized(task.minimizedConf);
				}
			});
		}
		
		public void waitForSpace() {
			tasks.waitForSpace();
		}
		
		public void waitForFinish() {
			tasks.waitForFinish();
		}
		
		public void cleanup() {
			
			// make sure all the tasks are finished before cleaning up
			// tasks that aren't finished yet won't cleanup properly and leak memory!
			tasks.waitForFinish();
			
			synchronized (pool) {
				
				// make sure everything has been returned to the pool
				if (pool.available() < pool.size()) {
					throw new Error(String.format("molecule pool in inconsistent state (only %d/%d molecules available), can't cleanup. this is a bug", pool.available(), pool.size()));
				}
				
				// cleanup the energy functions if needed
				for (EfuncAndPMC em : pool) {
					if (em.efunc instanceof EnergyFunction.NeedsCleanup) {
						((EnergyFunction.NeedsCleanup)em.efunc).cleanup();
					}
				}
				pool.clear();
			}
		}
	}
}
