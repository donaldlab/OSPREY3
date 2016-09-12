package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.AbstractTupleMatrix;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator.Result;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskListener;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.MoleculePool;
import edu.duke.cs.osprey.tools.Progress;

public class SimpleEnergyMatrixCalculator {
	
	private class MoleculeTask {
		
		public MoleculePool mols;
		
		protected Molecule getMolecule() {
			synchronized (mols) {
				return mols.checkout();
			}
		}
		
		protected void cleanup(Molecule mol, EnergyFunction efunc) {
			
			// release the molecule back to the pool
			synchronized (mols) {
				mols.release(mol);
			}
			
			// cleanup the energy function if needed
			if (efunc instanceof EnergyFunction.NeedsCleanup) {
				((EnergyFunction.NeedsCleanup)efunc).cleanup();
			}
		}
	}
	
	private class SingleTask extends MoleculeTask implements Runnable {
		
		public int pos1;
		public int rc1;
		public Result result;

		@Override
		public void run() {
			
			Molecule mol = getMolecule();
			EnergyFunction efunc = ecalc.getSingleEfunc(pos1, mol);
			RCTuple tup = new RCTuple();
			
			tup.set(pos1, rc1);
			result = ecalc.calc(efunc, tup, mol);
			
			cleanup(mol, efunc);
		}
	}
	
	private class PairTask extends MoleculeTask implements Runnable {
		
		public int pos1;
		public int rc1;
		public int pos2;
		public int numRcs2;
		public Result[] results;

		@Override
		public void run() {
			
			Molecule mol = getMolecule();
			EnergyFunction efunc = ecalc.getPairEfunc(pos1, pos2, mol);
			RCTuple tup = new RCTuple();
			
			results = new Result[numRcs2];
			for (int rc2=0; rc2<numRcs2; rc2++) {
				tup.set(pos1, rc1, pos2, rc2);
				results[rc2] = ecalc.calc(efunc, tup, mol);
			}
			
			cleanup(mol, efunc);
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
		EnergyMatrix emat = new EnergyMatrix(ecalc.getConfSpace(), Double.POSITIVE_INFINITY);
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
		
		// init molecule pool
		MoleculePool mols = new MoleculePool(ecalc.getConfSpace().m);
		
		// init task listeners
		TaskListener singleListener = new TaskListener() {
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
		};
		TaskListener pairListener = new TaskListener() {
			@Override
			public void onFinished(Runnable taskBase) {
				PairTask task = (PairTask)taskBase;
			
				for (int rc2=0; rc2<task.numRcs2; rc2++) {
					Result result = task.results[rc2];
					
					if (emat != null) {
						emat.setPairwise(task.pos1, task.rc1, task.pos2, rc2, result.getEnergy());
					}
					if (dofmat != null) {
						dofmat.setPairwise(task.pos1, task.rc1, task.pos2, rc2, result.getDofValues());
					}
				}
				
				progress.incrementProgress(task.results.length);
			}
		};
		
		System.out.println("Calculating energies with shell distribution: " + ecalc.getShellDistribution());
		
		for (int pos1=0; pos1<emat.getNumPos(); pos1++) {
			for (int rc1=0; rc1<emat.getNumConfAtPos(pos1); rc1++) {
			
				// NOTE: single terms tend to be much larger than pair terms,
				// so split up single terms into more different tasks than pair terms
				
				// singles
				SingleTask singleTask = new SingleTask();
				singleTask.mols = mols;
				singleTask.pos1 = pos1;
				singleTask.rc1 = rc1;
				tasks.submit(singleTask, singleListener);
				
				// pairs
				for (int pos2=0; pos2<pos1; pos2++) {
					
					PairTask pairTask = new PairTask();
					pairTask.mols = mols;
					pairTask.pos1 = pos1;
					pairTask.rc1 = rc1;
					pairTask.pos2 = pos2;
					pairTask.numRcs2 = emat.getNumConfAtPos(pos2);
					tasks.submit(pairTask, pairListener);
				}
			}
		}
		
		tasks.waitForFinish();
	}
}
