package edu.duke.cs.osprey.energy.compiled;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.gpu.Structs;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;

import java.util.List;
import java.util.Set;


/**
 * A basic conformation energy calculator implementation, written in pure Java.
 *
 * Not the fastest implementation available. Try the {@link NativeConfEnergyCalculator} first,
 * or the {@link CudaConfEnergyCalculator} if you have GPUs available.
 */
/*
 * But it is the simplest and easiest to understand implementation.
 * So if you're working on a new implemtation, this is a good place to start learning how things work.
 */
public class CPUConfEnergyCalculator implements ConfEnergyCalculator {

	public final ConfSpace confSpace;

	public CPUConfEnergyCalculator(ConfSpace confSpace) {
		this.confSpace = confSpace;
	}

	@Override
	public void close() {
		// nothing to do
	}

	@Override
	public ConfSpace confSpace() {
		return confSpace;
	}

	@Override
	public Structs.Precision precision() {
		return Structs.Precision.Float64;
	}

	@Override
	public EnergiedCoords calc(int[] conf, List<PosInter> inters) {

		// build the conformation coords
		AssignedCoords coords = confSpace.makeCoords(conf);

		double energy = calcEnergy(coords, inters);
		return new EnergiedCoords(coords, energy);
	}

	public double calcEnergy(AssignedCoords coords, List<PosInter> inters) {

		double energy = 0.0;

		for (PosInter inter : inters) {
			for (EnergyCalculator ecalc : confSpace.ecalcs) {
				energy += ecalc.calcEnergy(coords, inter);
			}
			energy += inter.weight*inter.offset;
		}

		return energy;
	}

	public double calcSubEnergy(AssignedCoords coords, List<PosInter> inters, Set<Integer> posIndices) {

		double energy = 0.0;

		for (PosInter inter : inters) {
			if (inter.isIncludedIn(posIndices)) {
				for (EnergyCalculator ecalc : confSpace.ecalcs) {
					energy += ecalc.calcEnergy(coords, inter);
				}
				energy += inter.weight*inter.offset;
			}
		}

		return energy;
	}

	@Override
	public EnergiedCoords minimize(int[] conf, List<PosInter> inters) {

		// build the conformation coords
		AssignedCoords coords = confSpace.makeCoords(conf);

		// TODO: can optimize by not including molecule rotation,translation DoFs
		//  unless interactions span across molecules

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
				return calcEnergy(coords, inters);
			}

			@Override
			public double getValForDOF(int dof, double val) {
				setDOF(dof, val);
				return calcSubEnergy(coords, inters, coords.dofs.get(dof).modifiedPosIndices());
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
