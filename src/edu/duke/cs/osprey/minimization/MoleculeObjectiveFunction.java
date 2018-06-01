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

package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleFactory1D;
import java.util.ArrayList;
import java.util.List;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.energy.EnergyFunction;

public class MoleculeObjectiveFunction implements ObjectiveFunction {

	private static final long serialVersionUID = -5301575611582359731L;

	public final ParametricMolecule pmol;
	public final EnergyFunction efunc;
	public final List<EnergyFunction> efuncsByDof;
	public final DoubleMatrix1D curDOFVals;

	public MoleculeObjectiveFunction(ParametricMolecule pmol, EnergyFunction efunc) {
		this.pmol = pmol;
		this.efunc = efunc;
		this.curDOFVals = DoubleFactory1D.dense.make(pmol.dofBounds.size());//will be set properly once we evaluate anything

		// init efunc if needed
		if (efunc instanceof EnergyFunction.NeedsInit) {
			((EnergyFunction.NeedsInit)efunc).init(pmol.mol, pmol.dofs, curDOFVals);
		}

		if (efunc instanceof EnergyFunction.DecomposableByDof) {
			efuncsByDof = ((EnergyFunction.DecomposableByDof)efunc).decomposeByDof(pmol.mol, pmol.dofs);
		} else {
			efuncsByDof = null;
		}
	}

	/**
	 * transition adapter, only here temporarily
	 */
	@Deprecated
	public MoleculeObjectiveFunction(MoleculeModifierAndScorer mof) {
		pmol = new ParametricMolecule(mof.getMolec(), mof.getDOFs(), new DofBounds(mof.getConstraints()));
		efunc = mof.getEfunc();
		efuncsByDof = new ArrayList<>();
		for (int d=0; d<pmol.dofs.size(); d++) {
			efuncsByDof.add(mof.getEfunc(d));
		}
		curDOFVals = mof.curDOFVals;//since the molecules are not deep-copied,
		//the DOF values shouldn't be either
	}

	public EnergyFunction getEfunc(int d) {
		if (efuncsByDof != null) {
			return efuncsByDof.get(d);
		}
		return efunc;
	}

	@Override
	public int getNumDOFs() {
		return pmol.dofs.size();
	}

	@Override
	public DoubleMatrix1D[] getConstraints() {
		return pmol.dofBounds.getBounds();
	}

	@Override
	public void setDOF(int d, double val) {
		curDOFVals.set(d, val);
		pmol.dofs.get(d).apply(val);
	}

	@Override
	public double getValForDOF(int d, double val) {
		setDOF(d, val);
		return getEfunc(d).getEnergy();
	}

	@Override
	public void setDOFs(DoubleMatrix1D x) {
		curDOFVals.assign(x);
		for (int d=0; d<x.size(); d++) {
			pmol.dofs.get(d).apply(x.get(d));
			//handleBlocksTogetherMaybe();//DEBUG!!!
		}
	}

	@Override
	public double getValue(DoubleMatrix1D x) {
		setDOFs(x);
		return efunc.getEnergy();
	}

	@Override
	public double getInitStepSize(int d) {
		return MoleculeModifierAndScorer.getInitStepSize(pmol.dofs.get(d));
	}

	@Override
	public boolean isDOFAngle(int d) {
		return MoleculeModifierAndScorer.isDOFAngle(pmol.dofs.get(d));
	}

	@Override
	public ArrayList<Integer> getInitFixableDOFs() {
		throw new UnsupportedOperationException("implement me!");
	}
}