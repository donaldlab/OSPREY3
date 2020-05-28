package edu.duke.cs.osprey.energy.compiled;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.parallelism.Parallelism;

import java.util.List;


public interface ConfEnergyCalculator extends AutoCloseable {

	class EnergiedCoords {

		public final AssignedCoords coords;
		public final double energy;
		/** non-null only when the conformation was minimized */
		public final DoubleMatrix1D dofValues;

		public EnergiedCoords(AssignedCoords coords, double energy, DoubleMatrix1D dofValues) {
			this.coords = coords;
			this.energy = energy;
			this.dofValues = dofValues;
		}

		public EnergiedCoords(AssignedCoords coords, double energy) {
			this(coords, energy, null);
		}
	}

	// don't make callers catch the Exception on close()
	@Override
	void close();


	ConfSpace confSpace();

	/**
	 * Build the conformation and calculate its rigid (ie unminimized) energy, using the provided interactions.
	 */
	EnergiedCoords calc(int[] conf, List<PosInter> inters);

	/**
	 * Minimize the conformation, using the provided interactions.
	 */
	EnergiedCoords minimize(int[] conf, List<PosInter> inters);

	default EnergiedCoords calcOrMinimize(int[] conf, List<PosInter> inters, boolean minimize) {
		if (minimize) {
			return minimize(conf, inters);
		} else {
			return calc(conf, inters);
		}
	}

	/**
	 * Calculate the rigid (ie unminimized) energy of the conformation, using the provided interactions.
	 */
	default double calcEnergy(int[] conf, List<PosInter> inters) {
		return calc(conf, inters).energy;
	}

	/**
	 * Calculate the minimized energy of the conformation, using the provided interactions.
	 */
	default double minimizeEnergy(int[] conf, List<PosInter> inters) {
		return minimize(conf, inters).energy;
	}

	default double calcOrMinimizeEnergy(int[] conf, List<PosInter> inters, boolean minimize) {
		if (minimize) {
			return minimizeEnergy(conf, inters);
		} else {
			return calcEnergy(conf, inters);
		}
	}

	default ConfSearch.EnergiedConf calcEnergy(ConfSearch.ScoredConf conf, List<PosInter> inters) {
		return new ConfSearch.EnergiedConf(conf, calcEnergy(conf.getAssignments(), inters));
	}

	default ConfSearch.EnergiedConf minimizeEnergy(ConfSearch.ScoredConf conf, List<PosInter> inters) {
		return new ConfSearch.EnergiedConf(conf, minimizeEnergy(conf.getAssignments(), inters));
	}

	default ConfSearch.EnergiedConf calcOrMinimizeEnergy(ConfSearch.ScoredConf conf, List<PosInter> inters, boolean minimize) {
		if (minimize) {
			return minimizeEnergy(conf, inters);
		} else {
			return calcEnergy(conf, inters);
		}
	}

	/**
	 * Builds the appropriate conformation energy calculator based on the desired parallelism.
	 */
	static ConfEnergyCalculator build(ConfSpace confSpace, Parallelism parallelism) {
		if (parallelism.numGpus > 0) {
			// TODO
			throw new UnsupportedOperationException("GPU energy calculation not implement yet");
		} else {
			return new CPUConfEnergyCalculator(confSpace);
		}
	}
}
