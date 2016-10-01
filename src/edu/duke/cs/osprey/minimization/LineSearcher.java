package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleMatrix1D;

public interface LineSearcher {
	
	void search(ObjectiveFunction f, DoubleMatrix1D x, int dof, DoubleMatrix1D mins, DoubleMatrix1D maxs);
}
