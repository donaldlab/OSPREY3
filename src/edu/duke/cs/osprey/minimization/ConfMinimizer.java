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
	
	private static class EfuncAndMol {
		public Molecule mol;
		public EnergyFunction efunc;
	}
	
	private class Task implements Runnable {
		
		public int index;
		public ScoredConf conf;
		public ConfSpace confSpace;
		public ObjectPool<EfuncAndMol> pool;
		public EnergiedConf minimizedConf;
		
		@Override
		public void run() {
			EfuncAndMol em = pool.checkout();
			minimizedConf = minimize(em.mol, conf, em.efunc, confSpace);
			pool.release(em);
		}
	}
	
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
	
	public List<EnergiedConf> minimize(List<ScoredConf> confs, Factory<? extends EnergyFunction,Molecule> efuncs, ConfSpace confSpace) {
		return minimize(confs, efuncs, confSpace, new TaskExecutor());
	}
	
	public List<EnergiedConf> minimize(List<ScoredConf> confs, Factory<? extends EnergyFunction,Molecule> efuncs, ConfSpace confSpace, TaskExecutor tasks) {
		
		// make space for all the minimized confs
		List<EnergiedConf> minimizedConfs = new ArrayList<>(confs.size());
		for (int i=0; i<confs.size(); i++) {
			minimizedConfs.add(null);
		}
		
		// init the listener to collect all the minimized confs
		Progress progress = new Progress(confs.size());
		TaskExecutor.TaskListener listener = new TaskExecutor.TaskListener() {
			@Override
			public void onFinished(Runnable taskBase) {
				
				Task task = (Task)taskBase;
				minimizedConfs.set(task.index, task.minimizedConf);
				progress.incrementProgress();
			}
		};
		
		// make a pool for molecules and energy functions
		// to keep concurrent tasks from racing each other
		ObjectPool<EfuncAndMol> pool = new ObjectPool<>(new Factory<EfuncAndMol,Void>() {
			@Override
			public EfuncAndMol make(Void context) {
				
				EfuncAndMol out = new EfuncAndMol();
				out.mol = new Molecule(confSpace.m);
				out.efunc = efuncs.make(out.mol);
				return out;
			}
		});
		
		// send all the minimization tasks
		for (int i=0; i<confs.size(); i++) {
			
			Task task = new Task();
			task.index = i;
			task.conf = confs.get(i);
			task.confSpace = confSpace;
			task.pool = pool;
			tasks.submit(task, listener);
		}
		
		tasks.waitForFinish();
		
		// cleanup the energy functions if needed
		for (EfuncAndMol em : pool) {
			if (em.efunc instanceof EnergyFunction.NeedsCleanup) {
				((EnergyFunction.NeedsCleanup)em.efunc).cleanup();
			}
		}
		
		return minimizedConfs;
	}
}
