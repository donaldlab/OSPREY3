package edu.duke.cs.osprey.energy.compiled;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;

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
	TaskExecutor tasks();

	/**
	 * Build the conformation and calculate its rigid (ie unminimized) energy, using the provided interactions.
	 */
	EnergiedCoords calc(int[] conf, List<PosInter> inters);

	/**
	 * Minimize the conformation, using the provided interactions.
	 */
	EnergiedCoords minimize(int[] conf, List<PosInter> inters);


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


	/**
	 * Builds the appropriate conformation energy calculator based on the desired parallelism.
	 */
	static ConfEnergyCalculator build(ConfSpace confSpace, Parallelism parallelism, TaskExecutor tasks) {
		if (parallelism.numGpus > 0) {
			// TODO
			throw new UnsupportedOperationException("GPU energy calculation not implement yet");
		} else {
			return new CPUConfEnergyCalculator(confSpace, tasks);
		}
	}
}
