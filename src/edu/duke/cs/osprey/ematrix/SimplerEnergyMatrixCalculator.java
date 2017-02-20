package edu.duke.cs.osprey.ematrix;

import java.io.File;

import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.control.Defaults;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.FFInterGen;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskListener;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Progress;

/**
 * Computes a matrix of energies between pairs of residue conformations to be used by A* search.
 */
public class SimplerEnergyMatrixCalculator {
	
	// NOTE: don't use GPUs on energy matrices, it's too slow
	// always use the CPU
	
	public static class Builder {
		
		/**
		 * The conformation space containing the strands to be designed.
		 * 
		 * If the strands are configured with continuous flexibility, the energy matrix will
		 * minimize residue conformation pairs before computing energies.
		 */
		private SimpleConfSpace confSpace;
		
		/** The forcefield parameters for energy calculation. */
		private ForcefieldParams ffparams;
		
		/** Available hardware for high-performance computation. */
		private Parallelism parallelism;
		private Factory<Minimizer,ObjectiveFunction> minimizerFactory;
		
		/**
		 * Path to file where energy matrix should be saved between computations.
		 * 
		 * @note Energy matrix computation can take a long time, but often the results
		 * can be reused between computations. Use a cache file to skip energy matrix
		 * computation on the next Osprey run if the energy matrix has already been
		 * computed once before.
		 * 
		 * @warning If design settings are changed between runs, Osprey will make
		 * some effort to detect that the energy matrix cache is out-of-date and compute a
		 * new energy matrix instead of usng the cached, incorrect one. Osprey might not detect
		 * all design changes though, and incorrectly reuse a cached energy matrix, so it
		 * is best to manually delete the entry matrix cache file after changing design settings.
		 */
		private File cacheFile = null;
		
		public Builder(SimpleConfSpace confSpace) {
			this.confSpace = confSpace;
			this.ffparams = Defaults.forcefieldParams;
			this.parallelism = Defaults.parallelism;
			this.minimizerFactory = (f) -> new SimpleCCDMinimizer(f);
		}
		
		public Builder setForcefieldParams(ForcefieldParams val) {
			ffparams = val;
			return this;
		}
		
		public Builder setParallelism(Parallelism val) {
			parallelism = val;
			return this;
		}
		
		public Builder setMinimizerFactory(Factory<Minimizer,ObjectiveFunction> val) {
			minimizerFactory = val;
			return this;
		}
		
		public SimplerEnergyMatrixCalculator build() {
			return new SimplerEnergyMatrixCalculator(confSpace, ffparams, parallelism.numThreads, minimizerFactory, cacheFile);
		}
	}
	
	public static Builder builder(SimpleConfSpace confSpace) {
		return new Builder(confSpace);
	}
	
	public static SimplerEnergyMatrixCalculator build(SimpleConfSpace confSpace) {
		return builder(confSpace).build();
	}
	
	private class SingleTask implements Runnable {
		
		public int pos1;
		public int rc1;
		public Minimizer.Result result;

		@Override
		public void run() {
			result = calcSingle(pos1, rc1);
		}
	}
	
	private class PairTask implements Runnable {
		
		public int pos1;
		public int rc1;
		public int pos2;
		public int numrc2;
		public Minimizer.Result[] results;

		@Override
		public void run() {
			results = new Minimizer.Result[numrc2];
			for (int rc2=0; rc2<numrc2; rc2++) {
				results[rc2] = calcPair(pos1, rc1, pos2, rc2);
			}
		}
	}
	
	private final SimpleConfSpace confSpace;
	private final int numThreads;
	private final Factory<Minimizer,ObjectiveFunction> minimizerFactory;
	private final File cacheFile;
	private final EnergyFunctionGenerator efuncgen;

	private SimplerEnergyMatrixCalculator(SimpleConfSpace confSpace, ForcefieldParams ffparams, int numThreads, Factory<Minimizer,ObjectiveFunction> minimizerFactory, File cacheFile) {
		this.confSpace = confSpace;
		this.numThreads = numThreads;
		this.minimizerFactory = minimizerFactory;
		this.cacheFile = cacheFile;
		
		efuncgen = new EnergyFunctionGenerator(ffparams);
	}
	
	public EnergyMatrix calcEnergyMatrix() {
		
		if (cacheFile != null) {
			return ObjectIO.readOrMake(
				cacheFile,
				EnergyMatrix.class,
				"energy matrix",
				(emat) -> emat.matches(confSpace),
				(context) -> calcEnergyMatrix()
			);
		} else {
			return reallyCalcEnergyMatrix();
		}
	}
	
	private EnergyMatrix reallyCalcEnergyMatrix() {
		
		// start the task executor
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(numThreads);
		
		// allocate the new matrix
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		
		// init task listeners
		Progress progress = new Progress(confSpace.getNumResConfs() + confSpace.getNumResConfPairs());
		TaskListener<SingleTask> singleListener = (taskBase) -> {
			SingleTask task = (SingleTask)taskBase;
			emat.setOneBody(task.pos1, task.rc1, task.result.energy);
			progress.incrementProgress();
		};
		TaskListener<PairTask> pairListener = (taskBase) -> {
			PairTask task = (PairTask)taskBase;
			for (int rc2=0; rc2<task.numrc2; rc2++) {
				emat.setPairwise(task.pos1, task.rc1, task.pos2, rc2, task.results[rc2].energy);
			}
			progress.incrementProgress(task.numrc2);
		};
		
		// send all the tasks
		System.out.println("Calculating energy matrix with " + progress.getTotalWork() + " entries...");
		for (int pos1=0; pos1<emat.getNumPos(); pos1++) {
			for (int rc1=0; rc1<emat.getNumConfAtPos(pos1); rc1++) {
			
				// NOTE: single terms tend to be much larger than pair terms,
				// so split up single terms into more different tasks than pair terms
				
				// singles
				SingleTask singleTask = new SingleTask();
				singleTask.pos1 = pos1;
				singleTask.rc1 = rc1;
				tasks.submit(singleTask, singleListener);
				
				// pairs
				for (int pos2=0; pos2<pos1; pos2++) {
					
					PairTask pairTask = new PairTask();
					pairTask.pos1 = pos1;
					pairTask.rc1 = rc1;
					pairTask.pos2 = pos2;
					pairTask.numrc2 = emat.getNumConfAtPos(pos2);
					tasks.submit(pairTask, pairListener);
				}
			}
		}
		
		tasks.waitForFinish();
		
		// cleanup
		tasks.cleanup();
		
		return emat;
	}
	
	public Minimizer.Result calcSingle(int pos1, int rc1) {
		RCTuple conf = new RCTuple(pos1, rc1);
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		EnergyFunction efunc = efuncgen.interactionEnergy(FFInterGen.makeIntraAndShell(confSpace, pos1, pmol.mol));
		return calcEnergy(pmol, conf, efunc);
	}
	
	public Minimizer.Result calcPair(int pos1, int rc1, int pos2, int rc2) {
		RCTuple conf = new RCTuple(pos1, rc1, pos2, rc2);
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		EnergyFunction efunc = efuncgen.interactionEnergy(FFInterGen.makeResPair(confSpace, pos1, pos2, pmol.mol));
		return calcEnergy(pmol, conf, efunc);
	}
	
	private Minimizer.Result calcEnergy(ParametricMolecule pmol, RCTuple conf, EnergyFunction efunc) {
		
		// no continuous DOFs? just evaluate the energy function
		if (pmol.dofs.isEmpty()) {
			return new Minimizer.Result(null, efunc.getEnergy());
		}
		
		// otherwise, minimize over the DOFs
		return minimizerFactory.make(new MoleculeObjectiveFunction(
			pmol,
			confSpace.makeBounds(conf),
			efunc
		)).minimize();
	}
}
