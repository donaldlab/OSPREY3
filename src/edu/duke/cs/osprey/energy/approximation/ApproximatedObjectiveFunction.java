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
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;

import java.util.ArrayList;
import java.util.List;


public class ApproximatedObjectiveFunction implements ObjectiveFunction {

	public interface Approximator {

		int numDofs();
		double getValue(DoubleMatrix1D x);
		double getValForDOF(int dof, double val, DoubleMatrix1D x);
		double error();

		interface Addable extends Approximator {
			List<Integer> dofBlockIds();
			List<Integer> dofCounts();
			int numParams();
			double train(List<Minimizer.Result> trainingSet, List<Minimizer.Result> testSet);
			void train(double value);
			Addable makeIdentity(List<Integer> dofBlockIds, List<Integer> dofCounts);
			void add(Addable other, double weight, double offset);
		}
	}


	public final Approximator approximator;
	public final ObjectiveFunction f;

	private final DoubleMatrix1D x;

	public ApproximatedObjectiveFunction(ObjectiveFunction f, Approximator approximator) {

		if (f.getNumDOFs() != approximator.numDofs()) {
			throw new IllegalArgumentException("number of degrees of freedom don't match!");
		}

		this.f = f;
		this.approximator = approximator;

		x = DoubleFactory1D.dense.make(f.getNumDOFs());
	}

	@Override
	public int getNumDOFs() {
		return f.getNumDOFs();
	}

	@Override
	public DoubleMatrix1D[] getConstraints() {
		return f.getConstraints();
	}

	@Override
	public void setDOFs(DoubleMatrix1D x) {
		this.x.assign(x);
		f.setDOFs(x);
	}

	@Override
	public void setDOF(int dof, double val) {
		this.x.set(dof, val);
		f.setDOF(dof, val);
	}

	@Override
	public double getValue(DoubleMatrix1D x) {
		this.x.assign(x);
		return f.getValue(x) + approximator.getValue(x);
	}

	@Override
	public double getValForDOF(int dof, double val) {
		this.x.set(dof, val);
		return f.getValForDOF(dof, val) + approximator.getValForDOF(dof, val, x);
	}

	@Override
	public double getInitStepSize(int dof) {
		return f.getInitStepSize(dof);
	}

	@Override
	public boolean isDOFAngle(int dof) {
		return f.isDOFAngle(dof);
	}

	@Override
	public ArrayList<Integer> getInitFixableDOFs() {
		return f.getInitFixableDOFs();
	}
}
