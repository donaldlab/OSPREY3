package edu.duke.cs.osprey.ematrix;

import java.io.File;

import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.FFInterGen;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Progress;

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
		private Parallelism parallelism = Parallelism.makeCpu(1);
		private Factory<Minimizer,ObjectiveFunction> minimizerFactory = (f) -> new SimpleCCDMinimizer(f);
		
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
		
		public Builder(SimpleConfSpace confSpace, ForcefieldParams ffparams) {
			this.confSpace = confSpace;
			this.ffparams = ffparams;
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
		
		public Builder setCacheFile(File val) {
			cacheFile = val;
			return this;
		}
		
		public SimplerEnergyMatrixCalculator build() {
			return new SimplerEnergyMatrixCalculator(confSpace, ffparams, parallelism.numThreads, minimizerFactory, cacheFile);
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
	
	/**
	 * Computes a matrix of energies between pairs of residue conformations to be used by A* search.
	 */
	public EnergyMatrix calcEnergyMatrix() {
		
		if (cacheFile != null) {
			return ObjectIO.readOrMake(
				cacheFile,
				EnergyMatrix.class,
				"energy matrix",
				(emat) -> emat.matches(confSpace),
				(context) -> reallyCalcEnergyMatrix()
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
		
		// send all the tasks
		Progress progress = new Progress(confSpace.getNumResConfs() + confSpace.getNumResConfPairs());
		System.out.println("Calculating energy matrix with " + progress.getTotalWork() + " entries...");
		for (int pos1=0; pos1<emat.getNumPos(); pos1++) {
			for (int rc1=0; rc1<emat.getNumConfAtPos(pos1); rc1++) {
			
				// NOTE: single terms tend to be much larger than pair terms,
				// so split up single terms into more different tasks than pair terms
				
				// singles
				final int fpos1 = pos1;
				final int frc1 = rc1;
				tasks.submit(
					() -> {
						return calcSingle(fpos1, frc1);
					},
					(Minimizer.Result result) -> {
						emat.setOneBody(fpos1, frc1, result.energy);
						progress.incrementProgress();
					}
				);
				
				// pairs
				for (int pos2=0; pos2<pos1; pos2++) {
					
					final int fpos2 = pos2;
					final int numrc2 = emat.getNumConfAtPos(pos2);
					tasks.submit(
						() -> {
							Minimizer.Result[] results = new Minimizer.Result[numrc2];
							for (int rc2=0; rc2<numrc2; rc2++) {
								results[rc2] = calcPair(fpos1, frc1, fpos2, rc2);
							}
							return results;
						},
						(Minimizer.Result[] results) -> {
							for (int rc2=0; rc2<numrc2; rc2++) {
								emat.setPairwise(fpos1, frc1, fpos2, rc2, results[rc2].energy);
							}
							progress.incrementProgress(numrc2);
						}
					);
				}
			}
		}
		
		tasks.waitForFinish();
		
		return emat;
	}
	
	/**
	 * Calculates a reference energy for each residue position and residue type
	 * based on the minimum energy of all residue conformations at that position
	 * and residue type.
	 */
	public SimpleReferenceEnergies calcReferenceEnergies() {
		
		// start the task executor
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		tasks.start(numThreads);
		
		SimpleReferenceEnergies eref = new SimpleReferenceEnergies();
		
		// send all the tasks
		Progress progress = new Progress(confSpace.getNumResConfs());
		System.out.println("Calculating reference energies for " + progress.getTotalWork() + " residue confs...");
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {
			
				int posi = pos.index;
				int rci = rc.index;
				String resType = rc.template.name;
				
				tasks.submit(
					() -> {
						return calcIntra(posi, rci);
					},
					(Minimizer.Result result) -> {
						
						// keep the min energy for each pos,resType
						Double e = eref.get(posi, resType);
						if (e == null || result.energy < e) {
							e = result.energy;
						}
						eref.set(posi, resType, e);
					}
				);
			}
		}
		
		tasks.waitForFinish();
		
		return eref;
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
	
	public Minimizer.Result calcIntra(int pos1, int rc1) {
		RCTuple conf = new RCTuple(pos1, rc1);
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		EnergyFunction efunc = efuncgen.interactionEnergy(FFInterGen.makeSingleRes(confSpace, pos1, pmol.mol));
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
