package edu.duke.cs.osprey.minimization;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.tools.Progress;

public class ConfMinimizer {
	
	public double minimize(Molecule mol, int[] conf, EnergyFunction efunc, ConfSpace confSpace) {
		RCTuple tuple = new RCTuple(conf);
		MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, confSpace, tuple, mol);
		Minimizer minimizer = new CCDMinimizer(mof, true);
		minimizer.minimize();
		mof.cleanup();
		return efunc.getEnergy();
	}
	
	public EnergiedConf minimize(Molecule mol, ScoredConf conf, EnergyFunction efunc, ConfSpace confSpace) {
		double energy = minimize(mol, conf.getAssignments(), efunc, confSpace);
		return new EnergiedConf(conf, energy);
	}
	
	public List<EnergiedConf> minimize(Molecule mol, List<ScoredConf> confs, EnergyFunction efunc, ConfSpace confSpace) {
		Progress progress = new Progress(confs.size());
		List<EnergiedConf> econfs = new ArrayList<>();
		for (ScoredConf conf : confs) {
			econfs.add(minimize(mol, conf, efunc, confSpace));
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
		
		Async async = new Async(efuncs, confSpace, tasks);
		async.setListener(new Async.Listener() {
			@Override
			public void onMinimized(EnergiedConf econf, Integer id) {
				econfs.set(id, econf);
				progress.incrementProgress();
			}
		});
		
		// minimize them all
		for (int i=0; i<confs.size(); i++) {
			async.minimizeAsync(confs.get(i), i);
		}
		async.waitForFinish();
		async.cleanup();
		
		return econfs;
	}
	
	public static class Async {
		
		private static class EfuncAndMol {
			public Molecule mol;
			public EnergyFunction efunc;
		}
		
		private class Task implements Runnable {
			
			public Integer id;
			public ScoredConf conf;
			public EnergiedConf minimizedConf;
			
			@Override
			public void run() {
				minimizedConf = minimizeSync(conf);
			}
		}
		
		public static interface Listener {
			void onMinimized(EnergiedConf econf, Integer id);
		}
	
		private ConfSpace confSpace;
		private TaskExecutor tasks;
		private Listener externalListener;
		private ConfMinimizer minimizer;
		private TaskExecutor.TaskListener taskListener;
		private ObjectPool<EfuncAndMol> pool;
	
		public Async(Factory<? extends EnergyFunction,Molecule> efuncs, ConfSpace confSpace, TaskExecutor tasks) {
			
			this.confSpace = confSpace;
			this.tasks = tasks;
			
			externalListener = null;
			minimizer = new ConfMinimizer();
			
			// init the listener to collect all the minimized confs
			taskListener = new TaskExecutor.TaskListener() {
				@Override
				public void onFinished(Runnable taskBase) {
					
					Task task = (Task)taskBase;
					
					// chain listeners
					externalListener.onMinimized(task.minimizedConf, task.id);
				}
			};
			
			// make a pool for molecules and energy functions
			// to keep concurrent tasks from racing each other
			pool = new ObjectPool<>(new Factory<EfuncAndMol,Void>() {
				@Override
				public EfuncAndMol make(Void context) {
					
					EfuncAndMol out = new EfuncAndMol();
					out.mol = new Molecule(confSpace.m);
					out.efunc = efuncs.make(out.mol);
					return out;
				}
			});
		}
		
		public EnergiedConf minimizeSync(ScoredConf conf) {
			EfuncAndMol em = pool.checkout();
			try {
				return minimizer.minimize(em.mol, conf, em.efunc, confSpace);
			} finally {
				pool.release(em);
			}
		}
		
		public void setListener(Listener val) {
			externalListener = val;
		}
		
		public void minimizeAsync(ScoredConf conf) {
			minimizeAsync(conf, null);
		}
		
		public void minimizeAsync(ScoredConf conf, Integer id) {
			
			if (externalListener == null) {
				throw new IllegalStateException("call setListener() before minimizng asynchronously, this is a bug");
			}
			
			// submit the minimization task
			Task task = new Task();
			task.id = id;
			task.conf = conf;
			tasks.submit(task, taskListener);
		}
		
		public void waitForFinish() {
			tasks.waitForFinish();
			
			// TEMP
			System.out.println("pool size: " + pool.size());
		}
		
		public void cleanup() {
			
			// cleanup the energy functions if needed
			for (EfuncAndMol em : pool) {
				if (em.efunc instanceof EnergyFunction.NeedsCleanup) {
					((EnergyFunction.NeedsCleanup)em.efunc).cleanup();
				}
			}
			pool.clear();
		}
	}
}
