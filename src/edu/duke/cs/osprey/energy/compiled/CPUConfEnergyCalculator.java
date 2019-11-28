package edu.duke.cs.osprey.energy.compiled;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;

import java.util.List;

import static edu.duke.cs.osprey.tools.Log.log;


public class CPUConfEnergyCalculator implements ConfEnergyCalculator {

	public final ConfSpace confSpace;

	public CPUConfEnergyCalculator(ConfSpace confSpace) {
		this.confSpace = confSpace;
	}

	@Override
	public ConfSpace confSpace() {
		return confSpace;
	}

	@Override
	public double calcEnergy(int[] conf, List<PosInter> inters) {

		// build the conformation coords
		AssignedCoords coords = confSpace.makeCoords(conf);

		double energy = 0.0;

		for (EnergyCalculator ecalc : confSpace.ecalcs) {
			for (PosInter inter : inters) {

				if (inter.posi1 == inter.posi2) {
					energy += inter.weight*ecalc.calcEnergySingle(coords, inter.posi1);
				} else if (inter.posi1 == PosInter.StaticPos) {
					energy += inter.weight*ecalc.calcEnergyStatic(coords, inter.posi2);
				} else if (inter.posi2 == PosInter.StaticPos) {
					energy += inter.weight*ecalc.calcEnergyStatic(coords, inter.posi1);
				} else {
					energy += inter.weight*ecalc.calcEnergyPair(coords, inter.posi1, inter.posi2);
				}
			}
		}

		return energy;
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
					energy += ecalc.calcEnergy(coords);
				}
				return energy;
			}

			@Override
			public double getValForDOF(int dof, double val) {

				setDOF(dof, val);

				double energy = 0.0;
				for (EnergyCalculator ecalc : confSpace.ecalcs) {
					energy += ecalc.calcSubEnergy(coords, coords.dofs.get(dof).modifiedPosIndices());
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

		return new EnergiedCoords(
			coords,
			result.energy,
			result.dofValues
		);
	}

	@Override
	public EnergiedCoords minimize(int[] conf) {
		return minimize(conf, PosInterDist.all(confSpace));
	}
}
