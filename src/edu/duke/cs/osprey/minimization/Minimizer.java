/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.minimization;

//Interface for minimizers.  Instantiated using an ObjectiveFunction

import cern.colt.matrix.DoubleMatrix1D;

/**
 *
 * @author mhall44
 */
public interface Minimizer {
	
    public static class Result {
    
        public DoubleMatrix1D dofValues;
        public double energy;
        
        public Result(DoubleMatrix1D dofValues, double energy) {
            this.dofValues = dofValues;
            this.energy = energy;
        }
    }
    
    Result minimize();
    
    public static interface NeedsCleanup extends Minimizer {
    	void cleanup();
    }
    
    public static interface Reusable extends Minimizer {
    	void init(ObjectiveFunction f);
    }
}
