package edu.duke.cs.osprey.energy.approximation;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.approximation.ApproximatedObjectiveFunction.Approximator;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Progress;

import java.io.File;
import java.util.*;

import static edu.duke.cs.osprey.tools.Log.log;


public class ApproximatorMatrixCalculator {

	public enum ApproximatorType {
		Quadratic
	}

	public final ConfEnergyCalculator confEcalc;

	/**
	 * Number of sample points per degree of freedom for the training set.
	 * Always use an odd number so we sample the center point of the voxel.
	 */
	private int numSamplesPerDoF = 9; // TODO: find a good default that balances speed and accuracy

	/**
	 * The type of model to use for approximating forcefield energies
	 */
	private ApproximatorType type = ApproximatorType.Quadratic;

	private File cacheFile = null;

	public ApproximatorMatrixCalculator(ConfEnergyCalculator confEcalc) {
		this.confEcalc = confEcalc;
	}

	public ApproximatorMatrixCalculator setNumSamplesPerDoF(int val) {
		numSamplesPerDoF = val;
		return this;
	}

	public ApproximatorMatrixCalculator setApproximatorType(ApproximatorType val) {
		type = val;
		return this;
	}

	public ApproximatorMatrixCalculator setCacheFile(File val) {
		cacheFile = val;
		return this;
	}

	public ApproximatorMatrix calc() {

		ApproximatorMatrix amat = new ApproximatorMatrix(confEcalc.confSpace);

		if (cacheFile != null && cacheFile.exists()) {
			amat.readFrom(cacheFile);
			log("read Approximator Matrix from file: %s", cacheFile.getAbsolutePath());
			return amat;
		}

		Progress progress = new Progress(confEcalc.confSpace.shellResNumbers.size()*confEcalc.confSpace.getNumResConfs());
		log("calculating %d approximators for %d RCs ...", progress.getTotalWork(), confEcalc.confSpace.getNumResConfs());

		for (String fixedResNum : confEcalc.confSpace.shellResNumbers) {
			for (SimpleConfSpace.Position pos1 : confEcalc.confSpace.positions) {
				for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {

					confEcalc.tasks.submit(
						() -> calc(fixedResNum, pos1.index, rc1.index),
						(approximator) -> {
							amat.set(pos1.index, rc1.index, fixedResNum, approximator);
							progress.incrementProgress();
						}
					);
				}
			}
		}

		confEcalc.tasks.waitForFinish();

		if (cacheFile != null) {
			amat.writeTo(cacheFile);
			log("wrote Approximator Matrix to file: %s", cacheFile.getAbsolutePath());
		}

		return amat;
	}

	private Approximator.Addable calc(String fixedResNum, int pos, int rc) {

		// make the residue interactions
		ResidueInteractions inters = new ResidueInteractions();
		inters.addPair(confEcalc.confSpace.positions.get(pos).resNum, fixedResNum);

		// make the molecule
		ParametricMolecule pmol = confEcalc.confSpace.makeMolecule(new RCTuple(pos, rc));
		int d = pmol.dofBounds.size();

		// make the energy function
		try (EnergyFunction ff = confEcalc.ecalc.makeEnergyFunction(pmol, inters)) {
			MoleculeObjectiveFunction f = new MoleculeObjectiveFunction(pmol, ff);

			DoubleMatrix1D x = DoubleFactory1D.dense.make(d);

			if (d > 0) {

				// have continuous flexibilty, sample dense from the config space

				int[] dims = new int[d];
				Arrays.fill(dims, numSamplesPerDoF);

				int numSamples = numSamplesPerDoF;
				for (int i=1; i<d; i++) {
					numSamples *= numSamplesPerDoF;
				}
				numSamples += 1;
				List<Minimizer.Result> samples = new ArrayList<>(numSamples);

				// sample points from a dense regular grid
				for (int[] p : new MathTools.GridIterable(dims)) {

					for (int d1=0; d1<d; d1++) {
						double min = pmol.dofBounds.getMin(d1);
						double max = pmol.dofBounds.getMax(d1);
						double xd = min + (max - min)*p[d1]/(dims[d1] - 1);
						x.set(d1, xd);
					}

					double energy = f.getValue(x);

					samples.add(new Minimizer.Result(x.copy(), energy));
				}

				// also sample the minimized center point
				samples.add(new SimpleCCDMinimizer(f).minimizeFromCenter());

				switch (type) {
					case Quadratic: return QuadraticApproximator.train(samples);
					default: throw new IllegalArgumentException("unknown approximator type: " + type);
				}

			} else {

				// no continuous flexibiltiy, just use the one energy value
				switch (type) {
					case Quadratic: return QuadraticApproximator.train(f.getValue(x));
					default: throw new IllegalArgumentException("unknown approximator type: " + type);
				}
			}
		}
	}
}
