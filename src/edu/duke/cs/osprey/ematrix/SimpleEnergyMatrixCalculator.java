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
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.ObjectPool.Checkout;
import edu.duke.cs.osprey.tools.Progress;

@Deprecated
public abstract class SimpleEnergyMatrixCalculator {
	
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
		
		System.out.println("Calculating energies...");
		
		for (int pos1=0; pos1<emat.getNumPos(); pos1++) {
			for (int rc1=0; rc1<emat.getNumConfAtPos(pos1); rc1++) {
			
				// NOTE: single terms tend to be much larger than pair terms,
				// so split up single terms into more different tasks than pair terms
				
				// singles
				final int fpos1 = pos1;
				final int frc1 = rc1;
				tasks.submit(
					() -> {
						try (Checkout<ParameterizedMoleculeCopy> pmol = pmols.autoCheckout()) {
							return ecalc.calcSingle(fpos1, frc1, pmol.get());
						}
					},
					(Minimizer.Result result) -> {
						if (emat != null) {
							emat.setOneBody(fpos1, frc1, result.energy);
						}
						if (dofmat != null) {
							dofmat.setOneBody(fpos1, frc1, result.dofValues);
						}
						progress.incrementProgress();
					}
				);
				
				// pairs
				for (int pos2=0; pos2<pos1; pos2++) {
					
					final int fpos2 = pos2;
					final int numrc2 = emat.getNumConfAtPos(pos2);
					tasks.submit(
						() -> {
							try (Checkout<ParameterizedMoleculeCopy> pmol = pmols.autoCheckout()) {
								Minimizer.Result[] results = new Minimizer.Result[numrc2];
								for (int rc2=0; rc2<numrc2; rc2++) {
									results[rc2] = ecalc.calcPair(fpos1, frc1, fpos2, rc2, pmol.get());
								}
								return results;
							}
						},
						(Minimizer.Result[] results) -> {
							if (emat != null) {
								for (int rc2=0; rc2<numrc2; rc2++) {
									emat.setPairwise(fpos1, frc1, fpos2, rc2, results[rc2].energy);
								}
							}
							if (dofmat != null) {
								for (int rc2=0; rc2<numrc2; rc2++) {
									dofmat.setPairwise(fpos1, frc1, fpos2, rc2, results[rc2].dofValues);
								}
							}
							progress.incrementProgress(numrc2);
						}
					);
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
