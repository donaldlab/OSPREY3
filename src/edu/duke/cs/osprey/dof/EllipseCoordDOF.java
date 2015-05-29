
package edu.duke.cs.osprey.dof;

import java.util.ArrayList;
import java.util.Arrays;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.colt.matrix.linalg.Algebra;

/**
 * @author am439
 *
 */
public class EllipseCoordDOF extends DegreeOfFreedom {

	// For a d-dimensional ellipse (i.e. d dihedral angles), the ellipse is
	// parametrized by a radius variable and d-1 dihedral angles. the radius
	// ranges from 0 to 1, the first d-2 dihedral angles range from 0 to pi,
	// and the d-1th dihedral angle ranges from 0 to 2pi. 
	
	boolean isRadius; 
	int index; // every angle is indexed
	double curVal; // current value of the parameter
	ArrayList<DegreeOfFreedom> dihedrals; // a vector of dihedrals in chi-space
									   // representing the residue being 
									   // parametrized. 
	DoubleMatrix2D A;
	DoubleMatrix1D c;
	double[] initVals; 
	boolean initSet = false; // have we set initial values yet

	public EllipseCoordDOF(
			boolean r,
			int a,
			double val,
			DoubleMatrix2D mat,
			DoubleMatrix1D center,
			ArrayList<DegreeOfFreedom> arrayList,
			double[] inits) {
		isRadius = r;
		index = a;
		curVal = val;
		A = mat;
		c = center;
		dihedrals = arrayList;
		initVals = inits;
	}
	
	@Override
	public void apply(double paramVal) {
		if (Math.abs(curVal-paramVal) < 0.001) { return; }
				
		if (!initSet) {
			for (int i=0; i<initVals.length; i++) { dihedrals.get(i).apply(initVals[i]); }
			initSet = true;
		}
				
		int dim = dihedrals.size(); // dimensionality
		double[] o = new double[dim]; //get old dihedrals
		for (int i=0; i<dim; i++) { o[i] = dihedrals.get(i).getCurVal(); }
//		System.out.print("Old dihedrals: ");
//		for (int i=0; i<o.length; i++) { System.out.print(o[i]+" "); }
//		System.out.println();

		double[] phi = getEllipsoidalCoords(o); // get old polars
		phi[index] = paramVal; //apply transformation
		double[] u = getCartesianCoords(phi);
		
		for (int i=0; i<u.length; i++) {
			if (Double.isNaN(u[i])) { u[i] = 0; }
		}
		
//		System.out.print("New dihedrals: ");
//		for (int i=0; i<u.length; i++) { System.out.print(u[i]+" "); }
//		System.out.println();
		
		for (int i=0; i<dihedrals.size(); i++) { dihedrals.get(i).apply(u[i]); }
		curVal = paramVal; // set value internally
	}

	public double[] getEllipsoidalCoords(double[] dihedrals) {
    	if (dihedrals.length==0) { return new double[0]; }
    	
		EigenvalueDecomposition evd = new EigenvalueDecomposition(A);
		DoubleMatrix2D Q = evd.getV();
		DoubleMatrix2D L = evd.getD();
		DoubleMatrix2D qT = Q.viewDice().copy();
		Algebra alg = new Algebra();
    	
    	// first transform the cartesian coordinates based on the ellipse
    	double[] s = new double[dihedrals.length];
    	for (int i=0; i<dihedrals.length; i++) { 
    		s[i] = dihedrals[i]-c.get(i);
    	}
    	double[] u = alg.mult(qT, DoubleFactory1D.dense.make(s)).toArray();
    	double[] x = new double[u.length];
    	for (int i=0; i<u.length; i++) {
    		x[i] = u[i]/Math.sqrt(L.get(i, i));
    	}    	
    	dihedrals = x;
    	
    	// now get elliptical coordinates
    	double radius = 0;
    	for (double d : dihedrals) { radius += d*d; }
    	radius = Math.sqrt(radius);
    	int n = dihedrals.length;
    	double[] phi = new double[n-1];
    	for (int i=0; i<n-1; i++) {
    		double d=0;
    		for (int j=i; j<n; j++) { d += dihedrals[j]*dihedrals[j]; }
    		double quot = dihedrals[i]/Math.sqrt(d);
    		phi[i] = Math.acos(quot);
    	}
    	if (dihedrals[n-1] < 0) { phi[n-2] = 2*Math.PI - phi[n-2]; }
    	double[] ellCoords = new double[n];
    	ellCoords[0] = radius;
    	for (int i=1; i<n; i++) { ellCoords[i] = phi[i-1]; }
    	return ellCoords;
    }
	
	public double[] getCartesianCoords(double[] polars) {
		if (polars.length==0) { return new double[0]; }
		
		EigenvalueDecomposition evd = new EigenvalueDecomposition(A);
		DoubleMatrix2D Q = evd.getV();
		DoubleMatrix2D L = evd.getD();
		DoubleMatrix2D qT = Q.viewDice().copy();
		Algebra alg = new Algebra();		
		
		int n = polars.length;		
		double radius = polars[0];
		double[] phis = new double[n-1];
		for (int i=1; i<n; i++) { phis[i-1] = polars[i]; }
		
		double[] cartCoords = new double[n];
		for (int i=0; i<n; i++) {
			double prod = 1;
			for (int j=0; j<i; j++) { prod = prod * Math.sin(phis[j]); }
			if (i<n-1) { prod = prod * Math.cos(phis[i]); } 			
			cartCoords[i] = radius*prod;
		}
		
		// transform cartesian coordinates back!
		
		double[] u = new double[cartCoords.length];
		for (int i=0; i<u.length; i++) {
			u[i] = cartCoords[i]*Math.sqrt(L.get(i, i));
		}
		double[] s = alg.mult(Q, DoubleFactory1D.dense.make(u)).toArray();
		double[] x = new double[s.length];
		for (int i=0; i<x.length; i++) { x[i] = s[i] + c.get(i); }
		
		return x;
	}

	public int getIndex() { return this.index; }
	
	public double getCurVal() { return this.curVal; }
	
}
