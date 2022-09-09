package edu.duke.cs.osprey.energy.compiled;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.gpu.Structs;
import edu.duke.cs.osprey.parallelism.Parallelism;

import java.util.List;


/**
 * The energy calculator interface for molecules created from compiled conformation spaces.
 *
 * Some implementations support multiple floating point precisions, depending on available hardware.
 */
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
	Structs.Precision precision();

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

	class MinimizationJob {

		public final int[] conf;
		public final List<PosInter> inters;

		public double energy;

		public MinimizationJob(int[] conf, List<PosInter> inters) {

			this.conf = conf;
			this.inters = inters;

			energy = Double.NaN;
		}
	}

	/**
	 * Calculate the minimized enegries of a batch of conformations and interactions.
	 * This will be the fastest minimization option by far on some implementations.
	 */
	default void minimizeEnergies(List<MinimizationJob> jobs) {
		for (var job : jobs) {
			job.energy = minimizeEnergy(job.conf, job.inters);
		}
	}

	default void calcEnergies(List<MinimizationJob> jobs) {
		for (var job : jobs) {
			job.energy = calcEnergy(job.conf, job.inters);
		}
	}

	default int maxBatchSize() {
		return 1;
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
	 * Builds the best conformation energy calculator based on the given resources.
	 */
	static ConfEnergyCalculator makeBest(ConfSpace confSpace, Parallelism parallelism, Structs.Precision precision) {

		// try GPUs first
		if (parallelism.numGpus > 0) {
			return new CudaConfEnergyCalculator(confSpace, precision, parallelism);
		}

		// TODO: prefer intel if the hardware suitably matched?

		// prefer the native ecalc, over the Java one
		return new NativeConfEnergyCalculator(confSpace, precision);
	}

	/**
	 * Builds the best conformation energy calculator based on the given resources,
	 * and using Float64 precision.
	 */
	static ConfEnergyCalculator makeBest(ConfSpace confSpace, Parallelism parallelism) {
		return makeBest(confSpace, parallelism, Structs.Precision.Float64);
	}
}
