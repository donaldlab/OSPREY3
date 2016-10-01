package edu.duke.cs.osprey.minimization;

import java.util.ArrayList;
import java.util.List;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;

public class SimpleCCDMinimizer implements Minimizer {
	
	private static final double MaxIterations = 30; // same as CCDMinimizer
	private static final double ConvergenceThreshold = 0.001; // same as CCDMinimizer
	
	private ObjectiveFunction ofunc;
	private List<Integer> dofs;
	private List<LineSearcher> lineSearchers;

	public SimpleCCDMinimizer(ObjectiveFunction ofunc) {
		this.ofunc = ofunc;
		this.lineSearchers = new ArrayList<>();
	}
	
	@Override
	public DoubleMatrix1D minimize() {
		
		// get bounds
		int numDofs = ofunc.getNumDOFs();
		DoubleMatrix1D mins = ofunc.getConstraints()[0];
		DoubleMatrix1D maxs = ofunc.getConstraints()[1];
		
		// set initial values to the center of the bounds
		DoubleMatrix1D vals = DoubleFactory1D.dense.make(numDofs);
		for (int i=0; i<numDofs; i++) {
			vals.set(i, (maxs.get(i) + mins.get(i))/2);
		}
		
		// which dofs should we update?
		dofs = new ArrayList<>();
		for (int i=0; i<numDofs; i++) {
			
			// is there even a range for this dof?
			if (mins.get(i) != maxs.get(i)) {
				dofs.add(i);
			}
		}
		
		// init the line searchers
		for (int i=0; i<dofs.size(); i++) {
			lineSearchers.add(new SurfingLineSearcher());
		}
		
		ccd(mins, maxs, vals);
		
		return vals;
	}
	
	private double ccd(DoubleMatrix1D mins, DoubleMatrix1D maxs, DoubleMatrix1D vals) {
		
		// ccd is pretty simple actually
		// just do a line search along each dimension until we stop improving
		// we deal with cycles by just capping the number of iterations
		
		// get the current objective function value
		double curf = ofunc.getValue(vals);
		
		for (int iter=0; iter<MaxIterations; iter++) {
			
			// update all the dofs using line search
			for (int i : dofs) {
				lineSearchers.get(i).search(ofunc, vals, i, mins, maxs);
			}
			
			// did we improve enough to keep going?
			double nextf = ofunc.getValue(vals);
			if (curf - nextf < ConvergenceThreshold) {
				
				// nope, we're done
				return nextf;
				
			} else {
				
				// yeah, keep going
				curf = nextf;
			}
		}
		
		return curf;
	}
}
