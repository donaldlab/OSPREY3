package edu.duke.cs.osprey.astar.conf.scoring;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;


public class HigherOrderGScorer implements AStarScorer {

	public final EnergyMatrix emat;

	public HigherOrderGScorer(EnergyMatrix emat) {
		this.emat = emat;
	}

	@Override
	public AStarScorer make() {
		return new HigherOrderGScorer(emat);
	}

	@Override
	public double calc(ConfIndex confIndex, RCs rcs) {

		// constant term
		// NOTE: Java's compiler is a bit dumb about closures,
		// so use a one-element array to work around it
		final double[] gscore = { emat.getConstTerm() };

		// one body energies
		for (int i=0; i<confIndex.numDefined; i++) {
			int pos1 = confIndex.definedPos[i];
			int rc1 = confIndex.definedRCs[i];

			gscore[0] += emat.getOneBody(pos1, rc1);
		}

		// pairwise energies
		for (int i=0; i<confIndex.numDefined; i++) {
			int pos1 = confIndex.definedPos[i];
			int rc1 = confIndex.definedRCs[i];

			for (int j=0; j<i; j++) {
				int pos2 = confIndex.definedPos[j];
				int rc2 = confIndex.definedRCs[j];

				gscore[0] += emat.getPairwise(pos1, rc1, pos2, rc2);
			}
		}

		// add higher order corrections
		int[] conf = Conf.make(confIndex);
		emat.forEachHigherOrderTupleIn(conf, (tuple, energy) -> {
			gscore[0] += energy;
		});

		return gscore[0];
	}
}
