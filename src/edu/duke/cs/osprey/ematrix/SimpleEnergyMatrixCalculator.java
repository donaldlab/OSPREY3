package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.AbstractTupleMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator.Result;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.MoleculePool;
import edu.duke.cs.osprey.tools.Progress;

public class SimpleEnergyMatrixCalculator {
	
	private class SingleTask implements Runnable {
		
		public MoleculePool mols;
		public int pos1;
		public int rc1;
		public Result result;

		@Override
		public void run() {
			Molecule mol = mols.checkout();
			result = ecalc.calcSingle(pos1, rc1, mol);
			mols.release(mol);
		}
	}
	
	private class PairTask implements Runnable {
		
		public MoleculePool mols;
		public int pos1;
		public int rc1;
		public int pos2;
		public int rc2;
		public Result result;

		@Override
		public void run() {
			Molecule mol = mols.checkout();
			result = ecalc.calcPair(pos1, rc1, pos2, rc2, mol);
			mols.release(mol);
		}
	}
	
	private abstract class TaskListener implements TaskExecutor.TaskListener {
		
		public EnergyMatrix emat;
		public DofMatrix dofmat;
		public Progress progress;
	}
	
	private class SingleListener extends TaskListener {
		
		@Override
		public void onFinished(Runnable taskBase) {
			SingleTask task = (SingleTask)taskBase;
			
			if (emat != null) {
				emat.setOneBody(task.pos1, task.rc1, task.result.getEnergy());
			}
			if (dofmat != null) {
				dofmat.setOneBody(task.pos1, task.rc1, task.result.getDofValues());
			}
			
			progress.incrementProgress();
		}
	}
	
	private class PairListener extends TaskListener {
		
		@Override
		public void onFinished(Runnable taskBase) {
			PairTask task = (PairTask)taskBase;
		
			if (emat != null) {
				emat.setPairwise(task.pos1, task.rc1, task.pos2, task.rc2, task.result.getEnergy());
			}
			if (dofmat != null) {
				dofmat.setPairwise(task.pos1, task.rc1, task.pos2, task.rc2, task.result.getDofValues());
			}
			
			progress.incrementProgress();
		}
	}
	
	private SimpleEnergyCalculator ecalc;
	
	public SimpleEnergyMatrixCalculator(SimpleEnergyCalculator ecalc) {
		this.ecalc = ecalc;
	}

	public EnergyMatrix calcEnergyMatrix() {
		return calcEnergyMatrix(null);
	}
	
	public EnergyMatrix calcEnergyMatrix(TaskExecutor tasks) {
		EnergyMatrix emat = new EnergyMatrix(ecalc.getConfSpace(), 0);
		calcMatrices(emat, null, tasks);
		return emat;
	}
	
	public DofMatrix calcDofMatrix() {
		return calcDofMatrix(null);
	}
	
	public DofMatrix calcDofMatrix(TaskExecutor tasks) {
		DofMatrix dofmat = new DofMatrix(ecalc.getConfSpace());
		calcMatrices(null, dofmat, tasks);
		return dofmat;
	}
	
	public void calcMatrices(EnergyMatrix emat, DofMatrix dofmat) {
		calcMatrices(emat, dofmat, null);
	}
	
	public void calcMatrices(EnergyMatrix emat, DofMatrix dofmat, TaskExecutor tasks) {
		
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
		
		// init task listeners
		SingleListener singleListener = new SingleListener();
		singleListener.emat = emat;
		singleListener.dofmat = dofmat;
		singleListener.progress = progress;
		
		PairListener pairListener = new PairListener();
		pairListener.emat = emat;
		pairListener.dofmat = dofmat;
		pairListener.progress = progress;
		
		// init molecule pool
		MoleculePool mols = new MoleculePool(ecalc.getConfSpace().m);
		
		System.out.println("Calculating energies with shell distribution: " + ecalc.getShellDistribution());
		for (int pos1=0; pos1<sizemat.getNumPos(); pos1++) {
			for (int rc1=0; rc1<sizemat.getNumConfAtPos(pos1); rc1++) {
				
				// singles
				SingleTask singleTask = new SingleTask();
				singleTask.mols = mols;
				singleTask.pos1 = pos1;
				singleTask.rc1 = rc1;
				tasks.submit(singleTask, singleListener);
				
				// pairs
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<sizemat.getNumConfAtPos(pos2); rc2++) {
					
						PairTask pairTask = new PairTask();
						pairTask.mols = mols;
						pairTask.pos1 = pos1;
						pairTask.rc1 = rc1;
						pairTask.pos2 = pos2;
						pairTask.rc2 = rc2;
						tasks.submit(pairTask, pairListener);
					}
				}
			}
		}
		
		tasks.waitForFinish();
	}
}
