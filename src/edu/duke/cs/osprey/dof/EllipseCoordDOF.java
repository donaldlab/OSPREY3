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

package edu.duke.cs.osprey.dof;

import java.util.ArrayList;
import java.util.Arrays;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.colt.matrix.linalg.Algebra;
import edu.duke.cs.osprey.tools.EllipseTransform;

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
	EllipseTransform et;

	public EllipseCoordDOF(
			boolean r,
			int a,
			double val,
			DoubleMatrix2D mat,
			ArrayList<DegreeOfFreedom> arrayList,
			double[] inits) {
		isRadius = r;
		index = a;
		curVal = val;
		A = mat;
		dihedrals = arrayList;
		initVals = inits;
	}

	@Override
	public void apply(double paramVal) {
		if (!initSet) { //set the center as the initial dihedrals
			for (int i=0; i<initVals.length; i++) { dihedrals.get(i).apply(initVals[i]); }
			double[] vals = new double[dihedrals.size()];
			for (int i=0; i<dihedrals.size(); i++) { vals[i] = dihedrals.get(i).getCurVal(); }
			c = DoubleFactory1D.dense.make(vals);
			initSet = true;
		}

		int dim = dihedrals.size(); // dimensionality
		double[] o = new double[dim]; //get old dihedrals
		for (int i=0; i<dim; i++) { o[i] = dihedrals.get(i).getCurVal(); }

		et = new EllipseTransform(A, c);
		double[] phi = et.getEllipsoidalCoords(o); // get old polars
		phi[index] = paramVal; //apply transformation
		double[] u = et.getCartesianCoords(phi);

		for (int i=0; i<u.length; i++) {
			if (Double.isNaN(u[i])) { u[i] = 0; }
		}

		for (int i=0; i<dihedrals.size(); i++) { dihedrals.get(i).apply(u[i]); }
		curVal = paramVal; // set value internally
	}

	public int getIndex() { return this.index; }

	public double getCurVal() { return this.curVal; }


	@Override
	public DOFBlock getBlock(){
		return null;
	}

	@Override
	public String getName() {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

}