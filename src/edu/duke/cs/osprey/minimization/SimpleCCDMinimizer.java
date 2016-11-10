package edu.duke.cs.osprey.minimization;

import java.util.ArrayList;
import java.util.List;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.tools.Factory;

public class SimpleCCDMinimizer implements Minimizer.NeedsCleanup {
	
	private static final double MaxIterations = 30; // same as CCDMinimizer
	private static final double ConvergenceThreshold = 0.001; // same as CCDMinimizer
	
	private ObjectiveFunction f;
	private List<ObjectiveFunction.OneDof> dofs;
	private List<LineSearcher> lineSearchers;

	public SimpleCCDMinimizer(ObjectiveFunction ofunc) {
		this(ofunc, new Factory<LineSearcher,Void>() {
			@Override
			public LineSearcher make(Void ignore) {
				return new SurfingLineSearcher();
			}
		});
	}
	
	public SimpleCCDMinimizer(ObjectiveFunction f, Factory<LineSearcher,Void> lineSearcherFactory) {
		
		this.f = f;
		
		// build the dofs
		dofs = new ArrayList<>();
		lineSearchers = new ArrayList<>();
		for (int d=0; d<f.getNumDOFs(); d++) {
			
			ObjectiveFunction.OneDof fd = new ObjectiveFunction.OneDof(f, d);
			dofs.add(fd);
			
			if (fd.getXMin() < fd.getXMax()) {
				LineSearcher lineSearcher = lineSearcherFactory.make(null);
				lineSearcher.init(fd);
				lineSearchers.add(lineSearcher);
			} else {
				lineSearchers.add(null);
			}
		}
	}
	
	@Override
	public DoubleMatrix1D minimize() {
		
		// init x to the center of the bounds
		int n = f.getNumDOFs();
		DoubleMatrix1D herex = DoubleFactory1D.dense.make(n);
		for (int d=0; d<n; d++) {
			ObjectiveFunction.OneDof dof = dofs.get(d);
			herex.set(d, (dof.getXMin() + dof.getXMax())/2);
		}
		
		DoubleMatrix1D nextx = herex.copy();
		
		// ccd is pretty simple actually
		// just do a line search along each dimension until we stop improving
		// we deal with cycles by just capping the number of iterations
		
		// get the current objective function value
		double herefx = f.getValue(herex);
		
		for (int iter=0; iter<MaxIterations; iter++) {
			
			// update all the dofs using line search
			for (int d=0; d<n; d++) {
				
				LineSearcher lineSearcher = lineSearchers.get(d);
				if (lineSearcher != null) {
					
					// get the next x value for this dof
					double xd = nextx.get(d);
					xd = lineSearcher.search(xd);
					nextx.set(d, xd);
				}
			}
			
			// how much did we improve?
			double nextfx = f.getValue(nextx);
			double improvement = herefx - nextfx;
			
			if (improvement > 0) {
				
				// take the step
				herex.assign(nextx);
				herefx = nextfx;
				
				if (improvement < ConvergenceThreshold) {
					break;
				}
				
			} else {
				break;
			}
		}
		
		return herex;
	}
	
	@Override
	public void cleanup() {
		for (LineSearcher lineSearcher : lineSearchers) {
			if (lineSearcher instanceof LineSearcher.NeedsCleanup) {
				((LineSearcher.NeedsCleanup)lineSearcher).cleanup();
			}
		}
	}
}
