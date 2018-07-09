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

import java.util.ArrayList;
import java.util.List;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.tools.Factory;

public class SimpleCCDMinimizer implements Minimizer.NeedsCleanup, Minimizer.Reusable {
	
	private static final double MaxIterations = 30; // same as CCDMinimizer
	private static final double ConvergenceThreshold = 0.001; // same as CCDMinimizer
	
	private Factory<LineSearcher,Void> lineSearcherFactory;
	private ObjectiveFunction f;
	private List<ObjectiveFunction.OneDof> dofs;
	private List<LineSearcher> lineSearchers;

	public SimpleCCDMinimizer() {
		this((context) -> new SurfingLineSearcher());
	}
	
	public SimpleCCDMinimizer(ObjectiveFunction f) {
		this();
		init(f);
	}
	
	public SimpleCCDMinimizer(Factory<LineSearcher,Void> lineSearcherFactory) {
		this.lineSearcherFactory = lineSearcherFactory;
		
		dofs = new ArrayList<>();
		lineSearchers = new ArrayList<>();
	}

	@Override
	public void init(ObjectiveFunction f) {
		
		this.f = f;
		
		// build the dofs
		dofs.clear();
		lineSearchers.clear();
		
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
	public Minimizer.Result minimizeFromCenter() {
		return minimizeFrom(f.getDOFsCenter());
	}

	@Override
	public Minimizer.Result minimizeFrom(DoubleMatrix1D startx) {

		int n = f.getNumDOFs();
		DoubleMatrix1D herex = startx.copy();
		DoubleMatrix1D nextx = startx.copy();

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

		// update the protein conf, one last time
		f.setDOFs(herex);

		return new Minimizer.Result(herex, herefx);
	}
	
	@Override
	public void clean() {
		for (LineSearcher lineSearcher : lineSearchers) {
			if (lineSearcher instanceof LineSearcher.NeedsCleanup) {
				((LineSearcher.NeedsCleanup)lineSearcher).cleanup();
			}
		}
	}
}
