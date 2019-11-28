package edu.duke.cs.osprey.energy.compiled;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;

import java.util.List;


public interface ConfEnergyCalculator {

	ConfSpace confSpace();

	/**
	 * Calculate the rigid (ie unminimized) energy of the conformation, using the provided interactions.
	 */
	double calcEnergy(int[] conf, List<PosInter> inters);

	/**
	 * Calculate the rigid (ie unminimized) energy of the conformation, using all interactions.
	 */
	default double calcEnergy(int[] conf) {
		return calcEnergy(conf, PosInterDist.all(confSpace()));
	}

	/**
	 * Calculate the minimized energy of the conformation, using the provided interactions.
	 */
	default double minimizeEnergy(int[] conf, List<PosInter> inters) {
		return minimize(conf, inters).energy;
	}

	/**
	 * Calculate the minimized energy of the conformation, using all interactions.
	 */
	default double minimizeEnergy(int[] conf) {
		return minimize(conf).energy;
	}

	class EnergiedCoords {

		public final AssignedCoords coords;
		public final double energy;
		public final DoubleMatrix1D dofValues;

		public EnergiedCoords(AssignedCoords coords, double energy, DoubleMatrix1D dofValues) {
			this.coords = coords;
			this.energy = energy;
			this.dofValues = dofValues;
		}
	}

	/**
	 * Minimize the conformation, using the provided interactions.
	 */
	EnergiedCoords minimize(int[] conf, List<PosInter> inters);

	/**
	 * Minimize the conformation, using using all interactions.
	 */
	EnergiedCoords minimize(int[] conf);
}
