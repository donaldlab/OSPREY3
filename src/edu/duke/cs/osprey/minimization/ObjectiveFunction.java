/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleMatrix1D;
import java.io.Serializable;

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
}

