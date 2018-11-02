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

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author mhall44
 */
//Objective function for CCDMinimizer or another general minimizer
public interface ObjectiveFunction extends Serializable {

    //Return number of degrees of freedom
    public int getNumDOFs();


    //Get constraints on the degrees of freedom
    public DoubleMatrix1D[] getConstraints();

	default DoubleMatrix1D getDOFsCenter() {
		int n = getNumDOFs();
		DoubleMatrix1D x = DoubleFactory1D.dense.make(n);
		for (int d=0; d<n; d++) {
			double xdmin = getConstraints()[0].get(d);
			double xdmax = getConstraints()[1].get(d);
			x.set(d, (xdmin + xdmax)/2);
		}
		return x;
	}

    //Set all DOFs to values in x (e.g. in the molecule)
    public void setDOFs(DoubleMatrix1D x);

    //Set just one degree of freedom
    public void setDOF(int dof, double val);

    //Value and gradient at a given point (specified as values for all DOFs)
    public double getValue(DoubleMatrix1D x);
    //public DoubleMatrix1D getGradient(DoubleMatrix1D x);

    //Value at a given value for a given DOF,
    //and, for efficiency, possibly omitting energy terms that don't depend on that DOF
    //Other DOFs kept as they are currently set
    //(omitted terms must be consistent for a given objective function and DOF)
    public double getValForDOF(int dof, double val);


    public double getInitStepSize(int dof);//Get a good initial step size for a DOF
    //(e.g. for first initial value checking in CCD)
    
    
    public boolean isDOFAngle(int dof);//Is the given degree of freedom an angle?
    //This is important because angles can wrap around (at 360-degree intervals)
    
    

    public ArrayList<Integer> getInitFixableDOFs();
    //If we're going to initialize full minimization with minimization over a limited number of DOFs,
    //these are the indices of the DOFs that will be fixed

    public static class OneDof implements Serializable {
		
		private ObjectiveFunction f;
		private int d;
		private double xdmin;
		private double xdmax;
		
		public OneDof(ObjectiveFunction f, int d) {
			this.f = f;
			this.d = d;
			this.xdmin = f.getConstraints()[0].get(d);
			this.xdmax = f.getConstraints()[1].get(d);
		}
		
		public ObjectiveFunction getParent() {
			return f;
		}
		
		public int getDimension() {
			return d;
		}
		
		public double getXMin() {
			return xdmin;
		}
		
		public double getXMax() {
			return xdmax;
		}
		
		public void setX(double xd) {
			f.setDOF(d, xd);
		}
		
		public double getValue(double xd) {
			return f.getValForDOF(d, xd);
		}
	}
	
	public static class DofBounds implements Serializable {
		
		private DoubleMatrix1D[] bounds;
		
		public static DofBounds concatenate(DofBounds ... bounds) {
			return concatenate(Arrays.asList(bounds));
		}
		
		public static DofBounds concatenate(List<DofBounds> subBoundsList) {
			
			// allocate the dof bounds
			int size = 0;
			for (DofBounds subBounds : subBoundsList) {
				if (subBounds != null) {
					size += subBounds.size();
				}
			}
			DofBounds bounds = new DofBounds(size);
			
			// copy the bounds
			int i=0;
			for (DofBounds subBounds : subBoundsList) {
				if (subBounds != null) {
					for (int j=0; j<subBounds.size(); j++) {
						bounds.set(i++, subBounds.getMin(j), subBounds.getMax(j));
					}
				}
			}
			return bounds;
		}
		
		public DofBounds(int numDofs) {
			bounds = new DoubleMatrix1D[] {
				DoubleFactory1D.dense.make(numDofs),
				DoubleFactory1D.dense.make(numDofs)
			};
		}
		
		public DofBounds(DofBounds other) {
			bounds = new DoubleMatrix1D[] {
				other.bounds[0].copy(),
				other.bounds[1].copy()
			};
		}
		
		public DofBounds(DoubleMatrix1D[] bounds) {
			this.bounds = bounds;
		}
		
		public DoubleMatrix1D[] getBounds() {
			return bounds;
		}
		
		public DoubleMatrix1D getMins() {
			return bounds[0];
		}
		
		public DoubleMatrix1D getMaxs() {
			return bounds[1];
		}
		
		public int size() {
			return getMins().size();
		}
	
		public double getMin(int d) {
			return getMins().get(d);
		}
		
		public double getMax(int d) {
			return getMaxs().get(d);
		}
		
		public double getCenter(int d) {
			return (getMin(d) + getMax(d))/2;
		}
		
		public void getCenter(DoubleMatrix1D out) {
			for (int d=0; d<size(); d++) {
				out.set(d, getCenter(d));
			}
		}

		public double getWidth(int d) {
			return getMax(d) - getMin(d);
		}
		
		public void set(int d, double min, double max) {
			bounds[0].set(d, min);
			bounds[1].set(d, max);
		}
		
		public double clamp(int d, double xd) {
			double min = getMin(d);
			double max = getMax(d);
			if (xd < min) {
				xd = min;
			} else if (xd > max) {
				xd = max;
			}
			return xd;
		}
		
		public void clamp(DoubleMatrix1D x) {
			for (int d=0; d<size(); d++) {
				x.set(d, clamp(d, x.get(d)));
			}
		}
		
		public double clampDelta(int d, double xd, double delta) {
			double min = getMin(d);
			double max = getMax(d);
			if (xd - delta < min) {
				delta = min - xd;
			} else if (xd + delta > max) {
				delta = max - xd;
			}
			return delta;
		}
		
		public void clampDelta(DoubleMatrix1D x, DoubleMatrix1D delta) {
			for (int d=0; d<size(); d++) {
				delta.set(d, clampDelta(d, x.get(d), delta.get(d)));
			}
		}
		
		public boolean isInBounds(int d, double xd) {
			return xd >= getMin(d) && xd <= getMax(d);
		}
		
		public boolean isInBounds(DoubleMatrix1D x) {
			for (int d=0; d<size(); d++) {
				if (!isInBounds(d, x.get(d))) {
					return false;
				}
			}
			return true;
		}

		@Override
		public String toString() {
			return toString(9, 4);
		}

		public String toString(int size, int precision) {
			String format = "%" + size + "." + precision + "f";
			format = "[" + format + "," + format + "]";
			StringBuilder buf = new StringBuilder();
			for (int d=0; d<size(); d++) {
				if (d > 0) {
					buf.append(", ");
				}
				buf.append(String.format(format, getMin(d), getMax(d)));
			}
			return buf.toString();
		}
	}

}
