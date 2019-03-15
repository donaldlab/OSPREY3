package edu.duke.cs.osprey.energy.approximation;


import cern.colt.matrix.DoubleMatrix1D;

public class NOPApproximator implements ApproximatedObjectiveFunction.Approximator {

	public final int d;

	public NOPApproximator(int d) {
		this.d = d;
	}

	@Override
	public int numDofs() {
		return d;
	}

	@Override
	public double getValue(DoubleMatrix1D x) {
		return 0;
	}

	@Override
	public double getValForDOF(int dof, double val, DoubleMatrix1D x) {
		return 0;
	}

	@Override
	public double error() {
		return 0;
	}
}
