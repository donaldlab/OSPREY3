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

package edu.duke.cs.osprey.tools;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import java.io.Serializable;

public class EllipseTransform implements Serializable {

	DoubleMatrix1D c;
	DoubleMatrix2D A;
	
	public EllipseTransform(DoubleMatrix2D mat, DoubleMatrix1D cent) {
		A = mat;
		c = cent;
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


}
