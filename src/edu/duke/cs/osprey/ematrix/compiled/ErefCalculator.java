package edu.duke.cs.osprey.ematrix.compiled;

import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.Progress;

import java.util.ArrayList;
import java.util.List;

import static edu.duke.cs.osprey.tools.Log.log;


/**
 * Calculate reference energies for a conformation space.
 *
 * This calculator uses just the 'single' energy of the design position as the reference energy.
 */
public class ErefCalculator {

	public static class Builder {

		public final ConfEnergyCalculator confEcalc;

		/**
		 * True to minimize conformations, false to use rigid conformations.
		 */
		private boolean minimize = true;

		public Builder(ConfEnergyCalculator confEcalc) {
			this.confEcalc = confEcalc;
		}

		public Builder setMinimize(boolean val) {
			minimize = val;
			return this;
		}

		public ErefCalculator build() {
			return new ErefCalculator(
				confEcalc,
				minimize
			);
		}
	}


	public final ConfEnergyCalculator confEcalc;
	public final boolean minimize;

	private ErefCalculator(ConfEnergyCalculator confEcalc, boolean minimize) {

		this.confEcalc = confEcalc;
		this.minimize = minimize;
	}

	public SimpleReferenceEnergies calc() {
		return calc(new TaskExecutor());
	}

	public SimpleReferenceEnergies calc(TaskExecutor tasks) {

		// allocate space
		SimpleReferenceEnergies eref = new SimpleReferenceEnergies();

		ConfSpace confSpace = confEcalc.confSpace();

		// count how much work there is to do
		Progress progress = new Progress(confSpace.countSingles());
		log("Calculating reference energies for %s position confs...", progress.getTotalWork());

		for (int posi=0; posi<confSpace.numPos(); posi++) {
			final int fposi = posi;
			for (int confi=0; confi<confSpace.numConf(posi); confi++) {
				final int fconfi = confi;
				tasks.submit(
					() -> {
						// use just the internal energy for the conformation
						List<PosInter> inters = new ArrayList<>();
						inters.add(new PosInter(fposi, fposi, 1.0, 0.0));

						int[] assignments = confSpace.assign(fposi, fconfi);
						return confEcalc.calcOrMinimizeEnergy(assignments, inters, minimize);
					},
					energy -> {

						// keep the min energy for each pos,resType
						String resType = confSpace.confType(fposi, fconfi);
						Double e = eref.get(fposi, resType);
						if (e == null || energy < e) {
							e = energy;
						}
						eref.set(fposi, resType, e);

						progress.incrementProgress();
					}
				);
			}
		}
		tasks.waitForFinish();

		return eref;
	}
}
