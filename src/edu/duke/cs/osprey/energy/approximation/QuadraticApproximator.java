package edu.duke.cs.osprey.energy.approximation;


import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.QRDecomposition;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.tools.IOable;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.List;


public class QuadraticApproximator implements ApproximatedObjectiveFunction.Approximator.Addable, IOable {

	public final int d;
	public final DoubleMatrix1D coefficients;

	private double maxe;

	public QuadraticApproximator(int d) {
		this.d = d;
		this.coefficients = DoubleFactory1D.dense.make(1 + d + d*(d + 1)/2);
		this.maxe = 0.0;
	}

	public QuadraticApproximator(int d, DoubleMatrix1D coefficients, double maxe) {
		this(d);
		this.coefficients.assign(coefficients);
		this.maxe = maxe;
	}

	@Override
	public int numDofs() {
		return d;
	}

	@Override
	public double getValue(cern.colt.matrix.DoubleMatrix1D x) {

		if (x.size() != d) {
			throw new IllegalArgumentException(String.format("x is wrong size (%d), expected %d", x.size(), d));
		}

		// constant term
		double v = coefficients.get(0);

		// linear terms
		for (int d1=0; d1<d; d1++) {
			v += coefficients.get(1 + d1)*x.get(d1);
		}

		// quadratic terms
		for (int d1=0; d1<d; d1++) {
			for (int d2=0; d2<=d1; d2++) {
				double c = coefficients.get(1 + d + d1*(d1 + 1)/2 + d2);
				v += c*x.get(d1)*x.get(d2);
			}
		}

		return v;
	}

	@Override
	public double getValForDOF(int dof, double val, DoubleMatrix1D x) {

		// constant term
		double v = coefficients.get(0);

		// linear terms
		for (int d1=0; d1<d; d1++) {
			double c = coefficients.get(1 + dof);
			double xd1 = d1 == dof ? val : x.get(d1);
			v += c*xd1;
		}

		// quadratic terms
		for (int d1=0; d1<d; d1++) {
			double xd1 = d1 == dof ? val : x.get(d1);
			for (int d2=0; d2<=d1; d2++) {
				double xd2 = d2 == dof ? val : x.get(d2);
				double c = coefficients.get(1 + d + d1*(d1 + 1)/2 + d2);
				v += c*xd1*xd2;
			}
		}

		return v;
	}

	@Override
	public double error() {
		return maxe;
	}

	@Override
	public QuadraticApproximator makeIdentity() {
		return new QuadraticApproximator(d);
	}

	@Override
	public void add(Addable other, double weight, double offset) {
		if (other instanceof QuadraticApproximator) {
			add((QuadraticApproximator)other, weight, offset);
		} else {
			throw new IllegalArgumentException("can't add different approximator types together:"
				+ "\n\t" + this.getClass().getName()
				+ "\n\t" + other.getClass().getName()
			);
		}
	}

	public void add(QuadraticApproximator other, double weight, double offset) {

		if (this.d != other.d) {
			throw new IllegalArgumentException(String.format("number of degrees of freedom don't match: %d != %d",
				this.d, other.d
			));
		}

		// add the weighted coefficients
		for (int i=0; i<this.coefficients.size(); i++) {
			this.coefficients.set(i, this.coefficients.get(i) + other.coefficients.get(i)*weight);
		}

		// add the offset
		this.coefficients.set(0, this.coefficients.get(0) + offset);

		this.maxe += other.maxe;
	}

	@Override
	public void writeTo(DataOutput out)
	throws IOException {
		for (int i=0; i<coefficients.size(); i++) {
			out.writeDouble(coefficients.get(i));
		}
		out.writeDouble(maxe);
	}

	@Override
	public void readFrom(DataInput in)
	throws IOException {
		for (int i=0; i<coefficients.size(); i++) {
			coefficients.set(i, in.readDouble());
		}
		maxe = in.readDouble();
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof QuadraticApproximator && equals((QuadraticApproximator)other);
	}

	public boolean equals(QuadraticApproximator other) {
		return this.d == other.d
			&& this.coefficients.equals(other.coefficients)
			&& this.maxe == other.maxe;
	}

	// TODO: use training and test sets?

	public static QuadraticApproximator train(List<Minimizer.Result> samples) {

		int d = samples.get(0).dofValues.size();

		// define the Ax=b linear system that describes our quadratic approximation
		DoubleMatrix2D b = DoubleFactory2D.dense.make(samples.size(), 1);
		DoubleMatrix2D A = DoubleFactory2D.dense.make(samples.size(), 1 + d + d*(d + 1)/2);

		for (int i=0; i<samples.size(); i++) {
			DoubleMatrix1D x = samples.get(i).dofValues;
			double energy = samples.get(i).energy;

			b.set(i, 0, energy);

			// degree 0 term
			A.set(i, 0, 1.0);

			// degree 1 terms
			for (int d1=0; d1<d; d1++) {
				A.set(i, 1 + d1, x.get(d1));
			}

			// degree 2 terms
			for (int d1=0; d1<d; d1++) {
				int d1offset = d1*(d1 + 1)/2;
				for (int d2=0; d2<=d1; d2++) {
					A.set(i, 1 + d + d1offset + d2, x.get(d1)*x.get(d2));
				}
			}
		}

		// solve Ax = b in least squares sense
		DoubleMatrix1D coefficients = new QRDecomposition(A).solve(b).viewColumn(0);

		// calculate the residual: Ax - b
		DoubleMatrix1D residual = new Algebra()
			.mult(A, coefficients)
			.assign(b.viewColumn(0), (ri, bi) -> ri - bi);

		// calculate the max training error?
		double maxe = 0.0;
		for (int i=0; i<samples.size(); i++) {
			maxe = Math.max(maxe, residual.get(i));
		}

		return new QuadraticApproximator(d, coefficients, maxe);
	}

	public static QuadraticApproximator train(double energy) {
		QuadraticApproximator approximator = new QuadraticApproximator(0);
		approximator.coefficients.set(0, energy);
		return approximator;
	}
}
