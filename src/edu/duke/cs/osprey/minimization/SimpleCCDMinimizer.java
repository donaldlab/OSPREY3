package edu.duke.cs.osprey.minimization;

import java.util.ArrayList;
import java.util.List;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.tools.Factory;

public class SimpleCCDMinimizer implements Minimizer {
	
	private static final double MaxIterations = 30; // same as CCDMinimizer
	private static final double ConvergenceThreshold = 0.001; // same as CCDMinimizer
	
	private ObjectiveFunction ofunc;
	private List<Integer> dofs;
	private List<LineSearcher> lineSearchers;
	private DoubleMatrix1D mins;
	private DoubleMatrix1D maxs;
	private DoubleMatrix1D vals;

	public SimpleCCDMinimizer(ObjectiveFunction ofunc) {
		this(ofunc, new Factory<LineSearcher,Void>() {
			@Override
			public LineSearcher make(Void context) {
				return new SurfingLineSearcher();
			}
		});
	}
	
	public SimpleCCDMinimizer(ObjectiveFunction ofunc, Factory<LineSearcher,Void> lineSearchers) {
		
		this.ofunc = ofunc;
		
		// get bounds
		int numDofs = ofunc.getNumDOFs();
		mins = ofunc.getConstraints()[0];
		maxs = ofunc.getConstraints()[1];
		
		// set initial values to the center of the bounds
		vals = DoubleFactory1D.dense.make(numDofs);
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
		this.lineSearchers = new ArrayList<>();
		for (int i=0; i<dofs.size(); i++) {
			this.lineSearchers.add(lineSearchers.make(null));
		}
	}
	
	@Override
	public DoubleMatrix1D minimize() {
		
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
				break;
				
			} else {
				
				// yeah, keep going
				curf = nextf;
			}
		}
		
		return vals;
	}
}
