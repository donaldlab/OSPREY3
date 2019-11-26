package edu.duke.cs.osprey.confspace.compiled;


import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.energy.compiled.EnergyCalculator;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;
import org.joml.Vector3d;
import org.joml.Vector3dc;

import java.util.ArrayList;
import java.util.List;


/**
 * A copy of the conf space atom coords with the desired assignments.
 */
public class AssignedCoords {

	public static int NotAssigned = -1;


	public final ConfSpace confSpace;

	/** conf indices for the design positions, in order */
	public final int[] assignments;

	private final int[] atomOffsetsByPos;

	/** atom coords for the conformations */
	private final CoordsList confCoords;

	/** degrees of freedom that modify the atom coords */
	public final List<DegreeOfFreedom> dofs;

	public AssignedCoords(ConfSpace confSpace, int[] assignments) {

		this.confSpace = confSpace;
		this.assignments = assignments;

		// copy all the conf coords into a single list

		// how many atoms do we need?
		atomOffsetsByPos = new int[confSpace.positions.length];
		int numCoords = 0;
		for (ConfSpace.Pos pos : confSpace.positions) {
			atomOffsetsByPos[pos.index] = numCoords;
			numCoords += pos.maxNumAtoms;
		}

		confCoords = new CoordsList(numCoords);
		dofs = new ArrayList<>();
		for (ConfSpace.Pos pos : confSpace.positions) {

			// get the conf, or skip this position if nothing was assigned
			int confi = assignments[pos.index];
			if (confi == NotAssigned) {
				continue;
			}

			// copy over the coords
			ConfSpace.Conf conf = pos.confs[confi];
			confCoords.copyFrom(conf.coords, atomOffsetsByPos[pos.index]);

			// make the motions and convert them into degrees of freedom
			for (int i=0; i<conf.motions.length; i++) {
				ContinuousMotion motion = conf.motions[i].build(this, pos);
				motion.appendDofs(dofs);
			}

			// TODO: add non-pos-specific DoFs?
		}
	}

	public double getStaticEnergy(int ffi) {
		return confSpace.staticEnergies[ffi];
	}

	public ConfSpace.IndicesSingle getIndices(int ffi, int posi) {

		// get the assignment, or null if nothing was assigned
		int confi = assignments[posi];
		if (confi == NotAssigned) {
			return null;
		}

		return confSpace.indicesSingles(ffi, posi, confi);
	}

	public ConfSpace.IndicesPair getIndices(int ffi, int posi1, int posi2) {

		// get the assignments, or null if nothing was assigned
		int confi1 = assignments[posi1];
		int confi2 = assignments[posi2];
		if (confi1 == NotAssigned || confi2 == NotAssigned) {
			return null;
		}

		return confSpace.indicesPairs(ffi, posi1, confi1, posi2, confi2);
	}

	public double[] getParams(int ffi, int paramsi) {
		return confSpace.ffparams(ffi, paramsi);
	}

	public void getStaticCoords(int atomi, Vector3d out) {
		confSpace.staticCoords.get(atomi, out);
	}

	public void getConfCoords(int posi, int atomi, Vector3d out) {
		confCoords.get(atomOffsetsByPos[posi] + atomi, out);
	}

	public void setConfCoords(int posi, int atomi, Vector3dc val) {
		confCoords.set(atomOffsetsByPos[posi] + atomi, val);
	}

	/**
	 * Calculate the energy of this conformation at the current
	 * atomic coordinates, ie without minimization.
	 */
	public double calcEnergy() {
		double energy = 0.0;
		for (EnergyCalculator ecalc : confSpace.ecalcs) {
			energy += ecalc.calcEnergy(this);
		}
		return energy;
	}

	/**
	 * Calculate the energy of this conformation at the current
	 * atomic coordinates (ie without minimization), only for the given positions.
	 */
	public double calcSubEnergy(int[] posIndices) {
		double energy = 0.0;
		for (EnergyCalculator ecalc : confSpace.ecalcs) {
			energy += ecalc.calcSubEnergy(this, posIndices);
		}
		return energy;
	}

	/**
	 * Calculate the energy of this conformation after minimizing
	 * from the starting conformation using the CCD minimization algorithm.
	 * After minimization, the atomic coordinates will be in the minimized positions.
	 */
	public Minimizer.Result minimize() {

		// create an objective function for minimization
		ObjectiveFunction f = new ObjectiveFunction() {

			@Override
			public int getNumDOFs() {
				return dofs.size();
			}

			@Override
			public DoubleMatrix1D[] getConstraints() {
				int n = dofs.size();
				DoubleMatrix1D mins = DoubleFactory1D.dense.make(n);
				DoubleMatrix1D maxs = DoubleFactory1D.dense.make(n);
				for (int d=0; d<n; d++) {
					mins.set(d, dofs.get(d).min());
					maxs.set(d, dofs.get(d).max());
				}
				return new DoubleMatrix1D[] { mins, maxs };
			}

			@Override
			public void setDOFs(DoubleMatrix1D x) {
				int n = dofs.size();
				for (int d=0; d<n; d++) {
					dofs.get(d).set(x.get(d));
				}
			}

			@Override
			public void setDOF(int dof, double val) {
				dofs.get(dof).set(val);
			}

			@Override
			public double getValue(DoubleMatrix1D x) {
				setDOFs(x);
				return calcEnergy();
			}

			@Override
			public double getValForDOF(int dof, double val) {
				setDOF(dof, val);
				return calcSubEnergy(dofs.get(dof).modifiedPosIndices());
			}

			@Override
			public double getInitStepSize(int dof) {
				return dofs.get(dof).initialStepSize();
			}
		};

		// minimize it!
		return new SimpleCCDMinimizer(f).minimizeFromCenter();
	}
}