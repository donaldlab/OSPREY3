/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleMatrix1D;
import java.io.Serializable;
import java.util.ArrayList;

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

    public static class OneDof {
		
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
	
	public static class DofBounds {
		
		private DoubleMatrix1D[] bounds;
		
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
	}

}

