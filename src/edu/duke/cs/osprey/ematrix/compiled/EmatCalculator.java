package edu.duke.cs.osprey.ematrix.compiled;

import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.tools.Progress;

import java.io.File;
import java.util.List;


/**
 * Calculates energy matrices for a conformation space.
 */
public class EmatCalculator {

	public static class Builder {

		public final ConfEnergyCalculator confEcalc;

		/**
		 * Defines how energies for single, pair, etc tuples should be distributed.
		 */
		private PosInterDist posInterDist = PosInterDist.DesmetEtAl1992;

		/**
		 * True to minimize conformations, false to use rigid conformations.
		 */
		private boolean minimize = true;

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

		public Builder setMinimize(boolean val) {
			minimize = val;
			return this;
		}

		public Builder setCacheFile(File val) {
			cacheFile = val;
			return this;
		}

		public EmatCalculator build() {
			return new EmatCalculator(
				confEcalc,
				posInterDist,
				minimize,
				cacheFile
			);
		}
	}


	public final ConfEnergyCalculator confEcalc;
	public final PosInterDist posInterDist;
	public final boolean minimize;
	public final File cacheFile;

	private EmatCalculator(ConfEnergyCalculator confEcalc, PosInterDist posInterDist, boolean minimize, File cacheFile) {

		this.confEcalc = confEcalc;
		this.posInterDist = posInterDist;
		this.minimize = minimize;
		this.cacheFile = cacheFile;
	}

	public EnergyMatrix calc() {

		// first, check the cache file for any previously-calculator energy matrices
		if (cacheFile != null) {
			// TODO: implement caching, key off all the arguments to this calculator
		}

		// not using cache, just calculate it
		return reallyCalc();

		// TODO: update the cache after calculating
	}

	private EnergyMatrix reallyCalc() {

		// allocate the new matrix
		EnergyMatrix emat = new EnergyMatrix(confEcalc.confSpace());

		ConfSpace confSpace = confEcalc.confSpace();

		// TODO: add reference energies?

		// count how much work there is to do
		// estimate work based on number of position interactions and the conf space size
		final int singleCost = posInterDist.single(confSpace, 0).size();
		final int pairCost = posInterDist.pair(confSpace, 0, 0).size();
		Progress progress = new Progress(
			confSpace.countConfSingles()*singleCost
				+ confSpace.countConfPairs()*pairCost
		);

		// TODO: use any parallelism provided by the confEcalc

		for (int posi1=0; posi1<emat.getNumPos(); posi1++) {
			for (int confi1=0; confi1<emat.getNumConfAtPos(posi1); confi1++) {

				// singles
				{
					int[] assignments = confSpace.assign(posi1, confi1);
					List<PosInter> inters = posInterDist.single(confSpace, posi1);
					double energy;
					if (minimize) {
						energy = confEcalc.minimizeEnergy(assignments, inters);
					} else {
						energy = confEcalc.calcEnergy(assignments, inters);
					}
					emat.setOneBody(posi1, confi1, energy);

					progress.incrementProgress(inters.size());
				}

				for (int posi2=0; posi2<posi1; posi2++) {
					for (int confi2=0; confi2<emat.getNumConfAtPos(posi2); confi2++) {

						// pairs
						int[] assignments = confSpace.assign(posi1, confi1, posi2, confi2);
						List<PosInter> inters = posInterDist.pair(confSpace, posi1, posi2);
						double energy;
						if (minimize) {
							energy = confEcalc.minimizeEnergy(assignments, inters);
						} else {
							energy = confEcalc.calcEnergy(assignments, inters);
						}
						emat.setPairwise(posi1, confi1, posi2, confi2, energy);

						progress.incrementProgress(inters.size());
					}
				}
			}
		}

		return emat;
	}
}
