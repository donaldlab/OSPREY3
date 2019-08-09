/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.energy.approximation;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.dof.DofInfo;
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
	 * Number of sample points per model parameter for the training set and test set
	 */
	private int numSamplesPerParam = 10;

	/**
	 * The type of model to use for approximating forcefield energies
	 */
	private ApproximatorType type = ApproximatorType.Quadratic;

	private File cacheFile = null;

	public ApproximatorMatrixCalculator(ConfEnergyCalculator confEcalc) {

		// build a new confEcalc based on the provided one,
		// but change the conf error budget to distibute the error budget among all the tuples evenly
		// so the errors for all the tuples can't add up in a single conf to exceed the budget
		int n = confEcalc.confSpace.positions.size();
		int maxTuplesPerConf = n + n*(n - 1)/2;
		this.confEcalc = new ConfEnergyCalculator.Builder(confEcalc.confSpace, confEcalc.ecalc)
			.setApproximationErrorBudget(confEcalc.approximationErrorBudget/maxTuplesPerConf)
			.setReferenceEnergies(confEcalc.eref)
			.setEnergyPartition(confEcalc.epart)
			.build();
	}

	public ApproximatorMatrixCalculator setNumSamplesPerParam(int val) {
		numSamplesPerParam = val;
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

		int numRCs = confEcalc.confSpace.getNumResConfs();
		Progress progress = new Progress(numRCs*(1 + confEcalc.confSpace.shellResNumbers.size()));
		log("calculating %d approximators for %d RCs ...", progress.getTotalWork(), numRCs);

		// singles and interactions with fixed residues
		for (SimpleConfSpace.Position pos1 : confEcalc.confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {

				// single
				confEcalc.tasks.submit(
					() -> calc(pos1, rc1),
					(approximator) -> {
						amat.set(pos1, rc1, approximator);
						progress.incrementProgress();
					}
				);

				// interactions with fixed residues
				for (String resNum : confEcalc.confSpace.shellResNumbers) {
					confEcalc.tasks.submit(
						() -> calc(pos1, rc1, resNum),
						(approximator) -> {
							amat.set(pos1, rc1, resNum, approximator);
							progress.incrementProgress();
						}
					);
				}
			}
		}

		// pairs
		for (SimpleConfSpace.Position pos1 : confEcalc.confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
				for (SimpleConfSpace.Position pos2 : confEcalc.confSpace.positions.subList(0, pos1.index)) {
					for (SimpleConfSpace.ResidueConf rc2: pos2.resConfs) {

						confEcalc.tasks.submit(
							() -> calc(pos1, rc1, pos2, rc2),
							(approximator) -> {
								amat.set(pos1, rc1, pos2, rc2, approximator);
								progress.incrementProgress();
							}
						);
					}
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

	public Approximator.Addable calc(SimpleConfSpace.Position pos, SimpleConfSpace.ResidueConf rc, String fixedResNum) {

		ResidueInteractions inters = new ResidueInteractions();
		inters.addPair(pos.resNum, fixedResNum);

		return calc(new RCTuple(pos.index, rc.index), inters);
	}

	public Approximator.Addable calc(SimpleConfSpace.Position pos, SimpleConfSpace.ResidueConf rc) {

		ResidueInteractions inters = new ResidueInteractions();
		inters.addPair(pos.resNum, pos.resNum);

		return calc(new RCTuple(pos.index, rc.index), inters);
	}

	public Approximator.Addable calc(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1, SimpleConfSpace.Position pos2, SimpleConfSpace.ResidueConf rc2) {

		ResidueInteractions inters = new ResidueInteractions();
		inters.addPair(pos1.resNum, pos2.resNum);

		return calc(new RCTuple(pos1.index, rc1.index, pos2.index, rc2.index), inters);
	}

	public Approximator.Addable calc(RCTuple tuple, ResidueInteractions inters) {

		// make the molecule
		ParametricMolecule pmol = confEcalc.confSpace.makeMolecule(tuple);

		// make the energy function
		try (EnergyFunction ff = confEcalc.ecalc.makeEnergyFunction(pmol, inters)) {
			MoleculeObjectiveFunction f = new MoleculeObjectiveFunction(pmol, ff);

			// gather the dof info
			DofInfo dofInfo = confEcalc.confSpace.makeDofInfo(tuple);

			// make the model
			Approximator.Addable approximator;
			switch (type) {
				case Quadratic: approximator = new QuadraticApproximator(dofInfo.ids, dofInfo.counts); break;
				default: throw new IllegalArgumentException("unknown approximator type: " + type);
			}

			if (pmol.dofs.isEmpty()) {

				// no continuous flexibiltiy, just use the one energy value
				approximator.train(ff.getEnergy());

			} else {

				// have continuous flexibilty, sample from the config space
				int numSamples = 1 + approximator.numParams()*numSamplesPerParam;
				Random rand = new Random(tuple.hashCode());
				approximator.train(
					sampleRandomly(pmol, f, numSamples, rand),
					sampleRandomly(pmol, f, numSamples, rand)
				);
			}

			return approximator;
		}
	}

	private List<Minimizer.Result> sampleRandomly(ParametricMolecule pmol, MoleculeObjectiveFunction f, int numSamples, Random rand) {

		List<Minimizer.Result> samples = new ArrayList<>(numSamples);

		// start with the minimized center point
		samples.add(new SimpleCCDMinimizer(f).minimizeFromCenter());

		// sample randomly
		for (int i=1; i<numSamples; i++) {

			DoubleMatrix1D x = DoubleFactory1D.dense.make(pmol.dofBounds.size());
			for (int d=0; d<pmol.dofBounds.size(); d++) {
				double min = pmol.dofBounds.getMin(d);
				double max = pmol.dofBounds.getMax(d);
				x.set(d, min + rand.nextDouble()*(max - min));
			}

			samples.add(new Minimizer.Result(x, f.getValue(x)));
		}

		// also sample the minimized center point
		samples.add(new SimpleCCDMinimizer(f).minimizeFromCenter());

		return samples;
	}

	private List<Minimizer.Result> sampleDensely(ParametricMolecule pmol, MoleculeObjectiveFunction f, int numSamplesPerDof) {

		int numDims = pmol.dofBounds.size();
		int[] dims = new int[numDims];
		Arrays.fill(dims, numSamplesPerDof);

		int numSamples = numSamplesPerDof;
		for (int i=1; i<numDims; i++) {
			numSamples *= numSamplesPerDof;
		}
		numSamples += 1;
		List<Minimizer.Result> samples = new ArrayList<>(numSamples);

		// start with the minimized center point
		samples.add(new SimpleCCDMinimizer(f).minimizeFromCenter());

		// sample points from a dense regular grid
		for (int[] p : new MathTools.GridIterable(dims)) {

			DoubleMatrix1D x = DoubleFactory1D.dense.make(numDims);
			for (int d=0; d<numDims; d++) {
				double min = pmol.dofBounds.getMin(d);
				double max = pmol.dofBounds.getMax(d);
				double xd = min + (max - min)*p[d]/(dims[d] - 1);
				x.set(d, xd);
			}

			samples.add(new Minimizer.Result(x, f.getValue(x)));
		}

		return samples;
	}
}
