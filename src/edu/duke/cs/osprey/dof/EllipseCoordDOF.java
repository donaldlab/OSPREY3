
package edu.duke.cs.osprey.dof;

import java.util.ArrayList;

import cern.colt.matrix.DoubleFactory1D;
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
	boolean isFinalTheta;
	double maxValue;
	double minValue = 0;
	double value; // current value of the parameter
	ArrayList<FreeDihedral> dihedrals; // a vector of dihedrals in chi-space
									   // representing the residue being 
									   // parametrized. 
	DoubleMatrix2D A;
	DoubleMatrix1D c;

	public EllipseCoordDOF(
			boolean r,
			boolean ft,
			int val,
			DoubleMatrix2D mat,
			DoubleMatrix1D center) {
		isRadius = r;
		isFinalTheta = ft;
		maxValue = (isRadius) ? 1 : (isFinalTheta) ? 2*Math.PI : Math.PI;
		value = val;
		A = mat;
		c = center;
	}
	
	/* (non-Javadoc)
	 * @see edu.duke.cs.osprey.dof.DegreeOfFreedom#apply(double)
	 */
	@Override
	public void apply(double paramVal) {
		int dim = dihedrals.size(); // dimensionality
		
		EigenvalueDecomposition evd = new EigenvalueDecomposition(A);
		DoubleMatrix2D Q = evd.getV();
		DoubleMatrix2D L = evd.getD();
		DoubleMatrix2D qT = Q.viewDice().copy();
		Algebra alg = new Algebra();
		

		// note: as of right now, either each dihedral will have to start
		// recording its value OR we need to add an arraylist of ellipsecorddofs 
		// in order to properly apply a given parameter. probably easier to just
		// add a tracking variable to each dihedral. this is because changing
		// one dihedral changes every ellipse coordinate, and changing one
		// ellipse coordinate changes every dihedral. 
		
		// get dihedrals
		double[] o = new double[dim];
		for (int i=0; i<dim; i++) { o[i] = dihedrals.get(i).getCurVal(); };
		// subtract center vector
		double[] x = new double[dim];
		for (int i=0; i<dim; i++) { x[i] = o[i] - c.get(i); }
		// rotate onto coordinate axes
		double[] s = alg.mult(qT, DoubleFactory1D.dense.make(x)).toArray();
		// scale to unit sphere
		double[] u = new double[dim];
		for (int i=0; i<dim; i++) { u[i] = s[i]/Q.get(i, i); }
		// transform into polar box
		double[] phi = new double[dim];
		// set radius
		for (int i=0; i<dim; i++) { phi[0] += u[i]*u[i]; }
		phi[0] = Math.sqrt(phi[0]);
		// set angles
		for (int i=0; i<dim; i++) { phi[i+1] = Math.acos(u[i]/phi[0]); }
		// divide the first d-2 angles by cos(value), the last one by sin(value)
		double[] a = phi.clone();
		for (int i=1; i<dim-1; i++) { a[i] = a[i]/Math.cos(value); }
		a[dim-1] = a[dim-1]/Math.sin(value);
		
		// now, multiply the first d-2 angles by cos(param< the last one by sin(param)
		for (int i=1; i<dim-1; i++) { a[i] = a[i]*Math.cos(paramVal); }
		a[dim-1] = a[dim-1] * Math.sin(paramVal);
		u = new double[dim];
		// define first n-1 cartesian vectors
		for (int i=0; i<dim-1; i++) {
			u[i] = phi[0]*Math.cos(phi[i+1]);
			for (int j=1; j<i+1; j++) { u[i] = u[i]*Math.sin(phi[j]); }
		}
		// define the nth cartesian vector
		u[dim-1] = phi[0]*Math.sin(phi[dim-1]);
		for (int j=1; j<dim-1; j++) { u[dim-1] = u[dim-1]*Math.sin(phi[j]); }
		// scale to ellipse 
		s = new double[dim];
		for (int i=0; i<dim; i++) { s[i] = u[i]*L.get(i, i); }
		// rotate the vector 
		double[] r = alg.mult(Q, DoubleFactory1D.dense.make(s)).toArray();
		// translate the vector
		double[] t = new double[dim];
		for (int i=0; i<dim; i++) { t[i] = r[i] + c.get(i); }
		
		// set all the dihedrals
		for (int i=0; i<dim; i++) { dihedrals.get(i).apply(t[i]); }
		
		// set the value
		value = paramVal;
	}

	public double getMaxValue() { return this.maxValue; }
	public double getMinValue() { return this.minValue; }
	
	
}
