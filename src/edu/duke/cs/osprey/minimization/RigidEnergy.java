package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleMatrix1D;

public class RigidEnergy extends HBFGSCCDMinimizer implements Minimizer {

	public RigidEnergy(ObjectiveFunction ofn, boolean useCorners) {
		super(ofn, useCorners);
	}

	public DoubleMatrix1D getDOFVals() {
        
        //First figure out the constraints and initial values
        //(since the minimizer might be used for many rotameric states, etc.,
        //this can't be done in the constructor)
        DoubleMatrix1D constr[] = objFcn.getConstraints();
        DOFmin = constr[0];
        DOFmax = constr[1];
        
        if(!compInitVals())//Initialize x
            return null;//No initial values found (this is currently only for IBEX)
        
        return x;
	}
	
	public DoubleMatrix1D minimize() {

		maxIter = 1;
		numIter = 1;
		
		return super.minimize();

	}

}
