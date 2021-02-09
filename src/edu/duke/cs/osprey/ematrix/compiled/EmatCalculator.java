package edu.duke.cs.osprey.ematrix.compiled;

import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.Progress;

import java.io.File;
import java.util.List;

import static edu.duke.cs.osprey.tools.Log.log;


/**
 * Calculates energy matrices for a conformation space.
 */
public class EmatCalculator {

	public static class Builder {

		public final ConfEnergyCalculator confEcalc;

		/** Defines what position interactions should be used for conformations and conformation fragments */
		private PosInterDist posInterDist = PosInterDist.DesmetEtAl1992;

		/** Reference energies used to design against the unfolded state, useful for GMEC designs */
		private SimpleReferenceEnergies eref = null;

		/**
		 * True to minimize conformations, false to use rigid conformations.
		 */
		private boolean minimize = true;

		/** True to include the static-static energies in conformation energies */
		private boolean includeStaticStatic = false;

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

		public Builder(ConfEnergyCalculator confEcalc) {
			this.confEcalc = confEcalc;
		}

		public Builder setPosInterDist(PosInterDist val) {
			posInterDist = val;
			return this;
		}

		public Builder setReferenceEnergies(SimpleReferenceEnergies val) {
			eref = val;
			return this;
		}

		public Builder setMinimize(boolean val) {
			minimize = val;
			return this;
		}

		public Builder setIncludeStaticStatic(boolean val) {
			includeStaticStatic = val;
			return this;
		}

		public Builder setCacheFile(File val) {
			cacheFile = val;
			return this;
		}

		public EmatCalculator build() {
			return new EmatCalculator(
				confEcalc,
				new PosInterGen(posInterDist, eref),
				minimize,
				includeStaticStatic,
				cacheFile
			);
		}
	}


	public final ConfEnergyCalculator confEcalc;
	public final PosInterGen posInterGen;
	public final boolean minimize;
	public final boolean includeStaticStatic;
	public final File cacheFile;

	private EmatCalculator(ConfEnergyCalculator confEcalc, PosInterGen posInterGen, boolean minimize, boolean includeStaticStatic, File cacheFile) {

		this.confEcalc = confEcalc;
		this.posInterGen = posInterGen;
		this.minimize = minimize;
		this.includeStaticStatic = includeStaticStatic;
		this.cacheFile = cacheFile;
	}

	public EnergyMatrix calc() {
		return calc(new TaskExecutor());
	}

	public EnergyMatrix calc(TaskExecutor tasks) {

		// first, check the cache file for any previously-calculated energy matrices
		if (cacheFile != null) {
			// TODO: implement caching, key off all the arguments to this calculator
		}

		// not using cache, just calculate it
		return reallyCalc(tasks);

		// TODO: update the cache after calculating
	}

	private EnergyMatrix reallyCalc(TaskExecutor tasks) {

		// allocate the new matrix
		EnergyMatrix emat = new EnergyMatrix(confEcalc.confSpace());

		ConfSpace confSpace = confEcalc.confSpace();

		// count how much work there is to do
		// estimate work based on number of position interactions and the conf space size
		final long staticCost = confSpace.avgAtomPairs(posInterGen.staticStatic());
		final long singleCost;
		final long pairCost;
		if (confSpace.numPos() > 0) {
			singleCost = confSpace.avgAtomPairs(posInterGen.single(confSpace, 0, 0));
			pairCost = confSpace.avgAtomPairs(posInterGen.pair(confSpace, 0, 0, 0, 0));
		} else {
			singleCost = 0;
			pairCost = 0;
		}
		int numSingles = confSpace.countSingles();
		int numPairs = confSpace.countPairs();
		Progress progress = new Progress(staticCost + numSingles*singleCost + numPairs*pairCost);
		log("Calculating energy matrix with %d entries", 1 + numSingles + numPairs);

		// static-static energy
		if (includeStaticStatic) {
			List<PosInter> inters = posInterGen.staticStatic();
			int[] conf = confSpace.assign();
			double energy = confEcalc.calcOrMinimizeEnergy(conf, inters, minimize);
			emat.setConstTerm(energy);
		}
		progress.incrementProgress(staticCost);

		for (int posi1=0; posi1<confSpace.numPos(); posi1++) {
			final int fposi1 = posi1;
			for (int confi1=0; confi1<confSpace.numConf(posi1); confi1++) {
				final int fconfi1 = confi1;

				// singles
				tasks.submit(
					() -> {
						int[] assignments = confSpace.assign(fposi1, fconfi1);
						List<PosInter> inters = posInterGen.single(confSpace, fposi1, fconfi1);
						return confEcalc.calcOrMinimizeEnergy(assignments, inters, minimize);
					},
					energy -> {
						emat.setOneBody(fposi1, fconfi1, energy);
						progress.incrementProgress(singleCost);
					}
				);

				for (int posi2=0; posi2<posi1; posi2++) {
					final int fposi2 = posi2;
					for (int confi2=0; confi2<confSpace.numConf(fposi2); confi2++) {
						final int fconfi2 = confi2;

						// pairs
						tasks.submit(
							() -> {
								int[] assignments = confSpace.assign(fposi1, fconfi1, fposi2, fconfi2);
								List<PosInter> inters = posInterGen.pair(confSpace, fposi1, fconfi1, fposi2, fconfi2);
								return confEcalc.calcOrMinimizeEnergy(assignments, inters, minimize);
							},
							energy -> {
								emat.setPairwise(fposi1, fconfi1, fposi2, fconfi2, energy);
								progress.incrementProgress(pairCost);
							}
						);
					}
				}
			}
		}
		tasks.waitForFinish();

		return emat;
	}
}
