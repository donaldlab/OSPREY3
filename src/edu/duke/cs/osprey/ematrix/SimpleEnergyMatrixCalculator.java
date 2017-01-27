package edu.duke.cs.osprey.ematrix;

import java.io.File;
import java.util.List;

import edu.duke.cs.osprey.confspace.AbstractTupleMatrix;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculePool;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskListener;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Progress;

@Deprecated
public abstract class SimpleEnergyMatrixCalculator {
	
	private class MoleculeTask {
		
		public ParameterizedMoleculePool pmols;
		
		protected ParameterizedMoleculeCopy checkout() {
			synchronized (pmols) {
				return pmols.checkout();
			}
		}
		
		protected void release(ParameterizedMoleculeCopy pmol) {
			synchronized (pmols) {
				pmols.release(pmol);
			}
		}
	}
	
	private class SingleTask extends MoleculeTask implements Runnable {
		
		public int pos1;
		public int rc1;
		public Minimizer.Result result;

		@Override
		public void run() {
			ParameterizedMoleculeCopy pmol = checkout();
			try {
				result = ecalc.calcSingle(pos1, rc1, pmol);
			} finally {
				release(pmol);
			}
		}
	}
	
	private class PairTask extends MoleculeTask implements Runnable {
		
		public int pos1;
		public int rc1;
		public int pos2;
		public int numrc2;
		public Minimizer.Result[] results;

		@Override
		public void run() {
			ParameterizedMoleculeCopy pmol = checkout();
			try {
				
				results = new Minimizer.Result[numrc2];
				for (int rc2=0; rc2<numrc2; rc2++) {
					results[rc2] = ecalc.calcPair(pos1, rc1, pos2, rc2, pmol);
				}
				
			} finally {
				release(pmol);
			}
		}
	}
	
	protected SimpleEnergyCalculator ecalc;
	protected TaskExecutor tasks;
	
	protected SimpleEnergyMatrixCalculator() {
		// only subclasses should directly make these
	}
	
	public EnergyMatrix calcEnergyMatrix(File cacheFile) {
		return ObjectIO.readOrMake(cacheFile, EnergyMatrix.class, "energy matrix", (context) -> calcEnergyMatrix());
	}
	
	public EnergyMatrix calcEnergyMatrix() {
		EnergyMatrix emat = new EnergyMatrix(ecalc.confSpace, Double.POSITIVE_INFINITY);
		calcMatrices(emat, null);
		return emat;
	}
	
	public DofMatrix calcDofMatrix() {
		DofMatrix dofmat = new DofMatrix(ecalc.confSpace);
		calcMatrices(null, dofmat);
		return dofmat;
	}
	
	public void calcMatrices(EnergyMatrix emat, DofMatrix dofmat) {
		
		if (emat != null && dofmat != null) {
			
			// make sure emat and dofmat match
			if (emat.getNumPos() != dofmat.getNumPos()) {
				throw new IllegalArgumentException("emat and dofmat must match size!");
			} else {
				for (int i=0; i<emat.getNumPos(); i++) {
					if (emat.getNumConfAtPos(i) != dofmat.getNumConfAtPos(i)) {
						throw new IllegalArgumentException("emat and dofmat must match size!");
					}
				}
			}
		}
		
		AbstractTupleMatrix<?> sizemat = null;
		if (emat != null) {
			sizemat = emat;
		}
		if (dofmat != null) {
			sizemat = dofmat;
		}
		if (sizemat == null) {
			throw new IllegalArgumentException("emat and dofmat cannot both be null");
		}
		
		// if we're not using fancy parallelism, then just do synchronous calculations
		if (tasks == null) {
			tasks = new TaskExecutor();
		}
		
		// count how much work there is to do
		long numWork = 0;
		for (int pos1=0; pos1<sizemat.getNumPos(); pos1++) {
			numWork += sizemat.getNumConfAtPos(pos1);
			for (int pos2=0; pos2<pos1; pos2++) {
				for (int rc1=0; rc1<sizemat.getNumConfAtPos(pos1); rc1++) {
					numWork += sizemat.getNumConfAtPos(pos2);
				}
			}
		}
		Progress progress = new Progress(numWork);
		
		// init molecule pool
		ParameterizedMoleculePool pmols = new ParameterizedMoleculePool(ecalc.confSpace);
		
		// init task listeners
		TaskListener singleListener = new TaskListener() {
			@Override
			public void onFinished(Runnable taskBase) {
				SingleTask task = (SingleTask)taskBase;
				
				if (emat != null) {
					emat.setOneBody(task.pos1, task.rc1, task.result.energy);
				}
				if (dofmat != null) {
					dofmat.setOneBody(task.pos1, task.rc1, task.result.dofValues);
				}
				
				progress.incrementProgress();
			}
		};
		TaskListener pairListener = new TaskListener() {
			@Override
			public void onFinished(Runnable taskBase) {
				PairTask task = (PairTask)taskBase;
			
				if (emat != null) {
					for (int rc2=0; rc2<task.numrc2; rc2++) {
						emat.setPairwise(task.pos1, task.rc1, task.pos2, rc2, task.results[rc2].energy);
					}
				}
				if (dofmat != null) {
					for (int rc2=0; rc2<task.numrc2; rc2++) {
						dofmat.setPairwise(task.pos1, task.rc1, task.pos2, rc2, task.results[rc2].dofValues);
					}
				}
				
				progress.incrementProgress(task.numrc2);
			}
		};
		
		System.out.println("Calculating energies...");
		
		for (int pos1=0; pos1<emat.getNumPos(); pos1++) {
			for (int rc1=0; rc1<emat.getNumConfAtPos(pos1); rc1++) {
			
				// NOTE: single terms tend to be much larger than pair terms,
				// so split up single terms into more different tasks than pair terms
				
				// singles
				SingleTask singleTask = new SingleTask();
				singleTask.pmols = pmols;
				singleTask.pos1 = pos1;
				singleTask.rc1 = rc1;
				tasks.submit(singleTask, singleListener);
				
				// pairs
				for (int pos2=0; pos2<pos1; pos2++) {
					
					PairTask pairTask = new PairTask();
					pairTask.pmols = pmols;
					pairTask.pos1 = pos1;
					pairTask.rc1 = rc1;
					pairTask.pos2 = pos2;
					pairTask.numrc2 = emat.getNumConfAtPos(pos2);
					tasks.submit(pairTask, pairListener);
				}
			}
		}
		
		tasks.waitForFinish();
	}
	
	public abstract void cleanup();
	
	
	public static class Cpu extends SimpleEnergyMatrixCalculator {
		
		private ThreadPoolTaskExecutor tasks;
		
		private Cpu(int numThreads, SimpleEnergyCalculator ecalc) {
			this.ecalc = ecalc;
			this.tasks = new ThreadPoolTaskExecutor();
			this.tasks.start(numThreads);
			super.tasks = tasks;
		}
		
		public Cpu(int numThreads, ForcefieldParams ffparams, ConfSpace confSpace, List<Residue> shellResidues) {
			this(numThreads, new SimpleEnergyCalculator.Cpu(ffparams, confSpace, shellResidues));
		}
		
		@Override
		public void cleanup() {
			tasks.stop();
		}
	}
	
	/**
	 * This is pretty slow compared to the CPU,
	 * so don't actually use it in the real world.
	 * Maybe someday it could be faster...
	 * Use the multi-threaded CPU calculator instead
	 */
	@Deprecated
	public static class Cuda extends SimpleEnergyMatrixCalculator {
		
		private GpuStreamPool pool;
		private ThreadPoolTaskExecutor tasks;
		
		public Cuda(int numGpus, int numStreamsPerGpu, ForcefieldParams ffparams, ConfSpace confSpace, List<Residue> shellResidues) {
			pool = new GpuStreamPool(numGpus, numStreamsPerGpu);
			ecalc = new SimpleEnergyCalculator.Cuda(pool, ffparams, confSpace, shellResidues);
			tasks = new ThreadPoolTaskExecutor();
			tasks.start(pool.getNumStreams());
			super.tasks = tasks;
		}
		
		@Override
		public void cleanup() {
			tasks.stop();
			pool.cleanup();
		}
	}
}
