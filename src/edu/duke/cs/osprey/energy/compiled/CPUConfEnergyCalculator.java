package edu.duke.cs.osprey.energy.compiled;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;

import java.util.List;


public class CPUConfEnergyCalculator implements ConfEnergyCalculator {

	public final ConfSpace confSpace;
	public final int numThreads;
	public final ThreadPoolTaskExecutor tasks;

	public CPUConfEnergyCalculator(ConfSpace confSpace, int numThreads) {

		this.confSpace = confSpace;
		this.numThreads = numThreads;

		tasks = new ThreadPoolTaskExecutor();
		tasks.start(numThreads);
	}

	@Override
	public void close() {
		tasks.close();
	}

	@Override
	public ThreadPoolTaskExecutor tasks() {
		return tasks;
	}

	@Override
	public ConfSpace confSpace() {
		return confSpace;
	}

	@Override
	public EnergiedCoords calc(int[] conf, List<PosInter> inters) {

		// build the conformation coords
		AssignedCoords coords = confSpace.makeCoords(conf);

		double energy = 0.0;

		for (EnergyCalculator ecalc : confSpace.ecalcs) {
			for (PosInter inter : inters) {
				energy += ecalc.calcEnergy(coords, inter);
			}
		}

		return new EnergiedCoords(coords, energy);
	}

	@Override
	public EnergiedCoords minimize(int[] conf, List<PosInter> inters) {

		// build the conformation coords
		AssignedCoords coords = confSpace.makeCoords(conf);

		// create an objective function for minimization
		ObjectiveFunction f = new ObjectiveFunction() {

			@Override
			public int getNumDOFs() {
				return coords.dofs.size();
			}

			@Override
			public DoubleMatrix1D[] getConstraints() {
				int n = coords.dofs.size();
				DoubleMatrix1D mins = DoubleFactory1D.dense.make(n);
				DoubleMatrix1D maxs = DoubleFactory1D.dense.make(n);
				for (int d=0; d<n; d++) {
					mins.set(d, coords.dofs.get(d).min());
					maxs.set(d, coords.dofs.get(d).max());
				}
				return new DoubleMatrix1D[] { mins, maxs };
			}

			@Override
			public void setDOFs(DoubleMatrix1D x) {
				int n = coords.dofs.size();
				for (int d=0; d<n; d++) {
					coords.dofs.get(d).set(x.get(d));
				}
			}

			@Override
			public void setDOF(int dof, double val) {
				coords.dofs.get(dof).set(val);
			}

			@Override
			public double getValue(DoubleMatrix1D x) {

				setDOFs(x);

				double energy = 0.0;
				for (EnergyCalculator ecalc : confSpace.ecalcs) {
					energy += ecalc.calcEnergy(coords, inters);
				}
				return energy;
			}

			@Override
			public double getValForDOF(int dof, double val) {

				setDOF(dof, val);

				double energy = 0.0;
				for (EnergyCalculator ecalc : confSpace.ecalcs) {
					energy += ecalc.calcSubEnergy(coords, inters, coords.dofs.get(dof).modifiedPosIndices());
				}
				return energy;
			}

			@Override
			public double getInitStepSize(int dof) {
				return coords.dofs.get(dof).initialStepSize();
			}
		};

		// minimize it!
		Minimizer.Result result = new SimpleCCDMinimizer(f).minimizeFromCenter();

		return new EnergiedCoords(coords, result.energy, result.dofValues);
	}
}
