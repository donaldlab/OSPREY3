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

	public final List<Integer> dofBlockIds;
	public final List<Integer> dofCounts;
	public final int numDofs;
	public final DoubleMatrix1D coefficients;

	private double maxe;

	private final int[] blockIndicesByDof;
	private final int[] dofOffsetsByBlock;

	public QuadraticApproximator(List<Integer> dofBlockIds, List<Integer> dofCounts) {

		this.dofBlockIds = dofBlockIds;
		this.dofCounts = dofCounts;
		this.numDofs = dofCounts.stream().mapToInt(i -> i).sum();
		this.coefficients = DoubleFactory1D.dense.make(1 + numDofs + numDofs*(numDofs + 1)/2);
		this.maxe = 0.0;

		// compute the dof<->block lookup tables
		blockIndicesByDof = new int[numDofs];
		dofOffsetsByBlock = new int[dofBlockIds.size()];
		int n = 0;
		for (int i=0; i<dofBlockIds.size(); i++) {
			dofOffsetsByBlock[i] = n;
			for (int j=0; j<dofCounts.get(i); j++) {
				blockIndicesByDof[n++] = i;
			}
		}
	}

	@Override
	public int numDofs() {
		return numDofs;
	}

	@Override
	public List<Integer> dofBlockIds() {
		return dofBlockIds;
	}

	@Override
	public List<Integer> dofCounts() {
		return dofCounts;
	}

	@Override
	public int numParams() {
		return coefficients.size();
	}

	@Override
	public double train(List<Minimizer.Result> trainingSet, List<Minimizer.Result> testSet) {

		// make sure all the samples have the right dimension
		for (Minimizer.Result sample : trainingSet) {
			if (sample.dofValues.size() != numDofs) {
				throw new IllegalArgumentException("samples have wrong number of dimensions");
			}
		}

		LinearSystem trainingSystem = new LinearSystem(trainingSet);
		LinearSystem testSystem = new LinearSystem(testSet);

		// solve Ax = b in least squares sense
		coefficients.assign(new QRDecomposition(trainingSystem.A).solve(trainingSystem.b).viewColumn(0));

		// calculate the residual (Ax - b) for the test set
		DoubleMatrix1D residual = new Algebra()
			.mult(testSystem.A, coefficients)
			.assign(testSystem.b.viewColumn(0), (ri, bi) -> Math.abs(ri - bi));

		// calculate the max error
		maxe = 0.0;
		for (int i=0; i<residual.size(); i++) {
			assert (residual.get(i) >= 0.0);
			maxe = Math.max(maxe, residual.get(i));
		}

		return maxe;
	}

	@Override
	public void train(double energy) {
		coefficients.set(0, energy);
		for (int i=1; i<coefficients.size(); i++) {
			coefficients.set(i, 0);
		}
	}

	private int index1(int d1) {
		return 1 + d1;
	}

	private int index2(int d1, int d2) {

		// make sure d2 <= d1
		if (d2 > d1) {
			int swap = d1;
			d1 = d2;
			d2 = swap;
		}

		return 1 + numDofs + d1*(d1 + 1)/2 + d2;
	}

	@Override
	public double getValue(DoubleMatrix1D x) {

		if (x.size() != numDofs) {
			throw new IllegalArgumentException(String.format("x is wrong size (%d), expected %d", x.size(), numDofs));
		}

		// constant term
		double v = coefficients.get(0);

		for (int d1=0; d1<numDofs; d1++) {

			// linear term
			double v1 = coefficients.get(index1(d1));

			// quadratic terms
			for (int d2=0; d2<=d1; d2++) {
				double c = coefficients.get(index2(d1, d2));
				v1 += c*x.get(d2);
			}

			v += v1*x.get(d1);
		}

		return v;
	}

	@Override
	public double getValForDOF(int d1, double val, DoubleMatrix1D x) {

		// linear term
		double v = coefficients.get(index1(d1));

		// quadratic terms
		for (int d2=0; d2<numDofs; d2++) {
			double x2 = x.get(d2);
			double c = coefficients.get(index2(d1, d2));
			v += c*x2;
		}

		v *= x.get(d1);

		// constant term
		v += coefficients.get(0);

		return v;
	}

	@Override
	public double error() {
		return maxe;
	}

	@Override
	public QuadraticApproximator makeIdentity(List<Integer> dofBlockIds, List<Integer> dofCounts) {
		return new QuadraticApproximator(dofBlockIds, dofCounts);
	}

	@Override
	public void add(Addable src, double weight, double offset) {
		if (src instanceof QuadraticApproximator) {
			add((QuadraticApproximator)src, this, weight, offset);
		} else {
			throw new IllegalArgumentException("can't add different approximator types together:"
				+ "\n\t" + this.getClass().getName()
				+ "\n\t" + src.getClass().getName()
			);
		}
	}

	public static void add(QuadraticApproximator src, QuadraticApproximator dst, double weight, double offset) {

		// match source blocks to destination blocks
		int[] dstBlockIndices = new int[src.dofBlockIds.size()];
		for (int srci=0; srci<src.dofBlockIds.size(); srci++) {
			int blockId = src.dofBlockIds.get(srci);

			int dsti = dst.dofBlockIds.indexOf(blockId);
			if (dsti < 0) {
				throw new IllegalArgumentException("destination approximator doesn't have dof block " + blockId);
			}

			// check the sizes just in case
			if (!dst.dofCounts.get(dsti).equals(src.dofCounts.get(srci))) {
				throw new IllegalArgumentException("block " + blockId + " has different sizes in different approximators");
			}

			dstBlockIndices[srci] = dsti;
		}

		// match source dofs to destination dofs
		int[] dstDofs = new int[src.numDofs];
		for (int srcd=0; srcd<src.numDofs; srcd++) {

			int srcb = src.blockIndicesByDof[srcd];

			int dofOffset = srcd - src.dofOffsetsByBlock[srcb];

			int dstb = dstBlockIndices[srcb];
			int dstd = dst.dofOffsetsByBlock[dstb] + dofOffset;

			dstDofs[srcd] = dstd;
		}

		// constant term (with the offset)
		dst.coefficients.set(0, dst.coefficients.get(0) + (src.coefficients.get(0) + offset)*weight);

		// linear terms
		for (int srcd1=0; srcd1<src.numDofs; srcd1++) {
			int dstd1 = dstDofs[srcd1];
			int srci = src.index1(srcd1);
			int dsti = dst.index1(dstd1);
			dst.coefficients.set(dsti, dst.coefficients.get(dsti) + src.coefficients.get(srci)*weight);
		}

		// quadratic terms
		for (int srcd1=0; srcd1<src.numDofs; srcd1++) {
			int dstd1 = dstDofs[srcd1];
			for (int srcd2=0; srcd2<=srcd1; srcd2++) {
				int dstd2 = dstDofs[srcd2];
				int srci = src.index2(srcd1, srcd2);
				int dsti = dst.index2(dstd1, dstd2);
				dst.coefficients.set(dsti, dst.coefficients.get(dsti) + src.coefficients.get(srci)*weight);
			}
		}

		// update the error
		dst.maxe += src.maxe;
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
		return this.dofBlockIds.equals(other.dofBlockIds)
			&& this.dofCounts.equals(other.dofCounts)
			&& this.coefficients.equals(other.coefficients)
			&& this.maxe == other.maxe;
	}

	// define the Ax=b linear system that describes our quadratic approximation
	private class LinearSystem {

		public final List<Minimizer.Result> samples;
		public final DoubleMatrix2D A;
		public final DoubleMatrix2D b;

		public LinearSystem(List<Minimizer.Result> samples) {

			this.samples = samples;

			A = DoubleFactory2D.dense.make(samples.size(), coefficients.size());
			b = DoubleFactory2D.dense.make(samples.size(), 1);

			for (int i=0; i<samples.size(); i++) {
				DoubleMatrix1D x = samples.get(i).dofValues;
				double energy = samples.get(i).energy;

				b.set(i, 0, energy);

				// degree 0 term
				A.set(i, 0, 1.0);

				// degree 1 terms
				for (int d1=0; d1<numDofs; d1++) {
					A.set(i, index1(d1), x.get(d1));
				}

				// degree 2 terms
				for (int d1=0; d1<numDofs; d1++) {
					for (int d2=0; d2<=d1; d2++) {
						A.set(i, index2(d1, d2), x.get(d1)*x.get(d2));
					}
				}
			}
		}
	}
}
