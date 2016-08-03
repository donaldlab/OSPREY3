package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleMatrix1D;

public class IntervalLineSearcher implements LineSearcher {
	
	private static final double Tau = (Math.sqrt(5) + 1)/2;
	private static final double MinIntervalWidth = 1e-3; // TODO: calibrate for different DOFs

	@Override
	public void search(ObjectiveFunction f, DoubleMatrix1D x, int dof, DoubleMatrix1D mins, DoubleMatrix1D maxs) {
		
		// use golden interval minimization
		// it won't necessarily pick the closest minimum,
		// but it will quickly converge to one of them arbitrarily
		// https://en.wikipedia.org/wiki/Golden_section_search
		
		// start with the interval equal to the full range
		double a = mins.get(dof);
		double b = maxs.get(dof);
		
		while (Math.abs(a - b) > MinIntervalWidth) {
			
			// use the golden ratio to sample points in the interior of the interval
			double step = (b - a)/Tau;
			double c = b - step;
			double d = a + step;
			
			// update interval bounds based on which side is lower
			double fc = f.getValForDOF(dof, c);
			double fd = f.getValForDOF(dof, d);
			if (fc < fd) {
				b = d;
			} else {
				a = c;
			}
		}
		
		// the interval is small now, pick the center as the final value
		x.set(dof, (a + b)/2);
	}

	@Override
	public void search(ObjectiveFunction f, DoubleMatrix1D x, DoubleMatrix1D vec, DoubleMatrix1D mins, DoubleMatrix1D maxs) {
		// TODO: line search in arbitrary direction
		throw new Error("implement me");
	}
}
