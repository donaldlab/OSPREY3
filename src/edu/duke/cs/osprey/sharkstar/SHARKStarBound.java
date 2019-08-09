package edu.duke.cs.osprey.sharkstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.parallelism.Parallelism;

import java.math.BigDecimal;

public class SHARKStarBound extends MARKStarBound {

	/**
	 * Constructor to make a default MARKStarBound Class
	 * @param confSpace the partition function conformation space
	 * @param rigidEmat the rigid pairwise energy matrix
	 * @param minimizingEmat the parwise-minimized energy matrix
	 * @param minimizingConfEcalc the energy calculator to calculate minimized conf energies
	 * @param rcs information on possible rotamers at all design positions
	 * @param parallelism information for threading
	 */
	public SHARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
						 ConfEnergyCalculator minimizingConfEcalc, RCs rcs, Parallelism parallelism) {
	    super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);
    }

	public SHARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
						  ConfEnergyCalculator minimizingConfEcalc, RCs rcs, Parallelism parallelism,
						  SHARKStarBound precomputedFlex) {
		super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);

		/*
		Now we need to do a couple things.
		For now, let's assume we are working with a single sequence

		TODO: Pass the precomputed tree to the new class
		 */

		super.rootNode = precomputedFlex.rootNode;

		/*
		TODO: Go through the tree to make sure that the assignments are compatible with the new confspace

		TODO: Make sure all of the energy matrices (including the correction emats) are compatible with the new confspace

		TODO: Do something with score / object contexts?

		TODO: Go through and update bounds

		TODO: Populate queue
		 */
	}

	public BigDecimal getLowerBound (Sequence seq){
		throw new UnsupportedOperationException("getLowerBound(seq) is not yet implemented");
	}

	public BigDecimal getUpperBound (Sequence seq){
		throw new UnsupportedOperationException("getUpperBound(seq) is not yet implemented");
	}

}
