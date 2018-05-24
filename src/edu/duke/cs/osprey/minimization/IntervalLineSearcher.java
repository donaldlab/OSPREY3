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

public class IntervalLineSearcher implements LineSearcher {
	
	private static final double Tau = (Math.sqrt(5) + 1)/2;
	private static final double MinIntervalWidth = 1e-3; // TODO: calibrate for different DOFs
	
	private ObjectiveFunction.OneDof f;

	@Override
	public void init(ObjectiveFunction.OneDof f) {
		this.f = f;
	}

	@Override
	public double search(double xd) {
		
		// use golden interval minimization
		// it won't necessarily pick the closest minimum,
		// but it will converge to one of them arbitrarily
		// https://en.wikipedia.org/wiki/Golden_section_search
		
		// start with the interval equal to the full range
		double a = f.getXMin();
		double b = f.getXMax();
		
		while (Math.abs(a - b) > MinIntervalWidth) {
			
			// use the golden ratio to sample points in the interior of the interval
			double step = (b - a)/Tau;
			double c = b - step;
			double d = a + step;
			
			// update interval bounds based on which side is lower
			double fc = f.getValue(d);
			double fd = f.getValue(d);
			if (fc < fd) {
				b = d;
			} else {
				a = c;
			}
		}
		
		// the interval is small now, pick the center as the final value
		double x = (a + b)/2;
		f.setX(x);
		return x;
	}
}
