package edu.duke.cs.osprey.ematrix.compiled;

import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
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

		// allocate space
		SimpleReferenceEnergies eref = new SimpleReferenceEnergies();

		ConfSpace confSpace = confEcalc.confSpace();

		// count how much work there is to do
		Progress progress = new Progress(confSpace.countSingles());
		log("Calculating reference energies for %s position confs...", progress.getTotalWork());

		// TODO: use any parallelism provided by the confEcalc

		for (int posi=0; posi<confSpace.numPos(); posi++) {
			for (int confi=0; confi<confSpace.numConf(posi); confi++) {

				// use just the internal energy for the conformation
				List<PosInter> inters = new ArrayList<>();
				inters.add(new PosInter(posi, posi, 1.0, 0.0));

				int[] assignments = confSpace.assign(posi, confi);
				double energy = confEcalc.calcOrMinimizeEnergy(assignments, inters, minimize);

				// keep the min energy for each pos,resType
				String resType = confSpace.confType(posi, confi);
				Double e = eref.get(posi, resType);
				if (e == null || energy < e) {
					e = energy;
				}
				eref.set(posi, resType, e);

				progress.incrementProgress();
			}
		}

		return eref;
	}
}
