/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.minimization;

//Interface for minimizers.  Instantiated using an ObjectiveFunction

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.tools.AutoCleanable;

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
    
    public static interface NeedsCleanup extends Minimizer, AutoCleanable {}
    
    public static interface Reusable extends Minimizer {
    	void init(ObjectiveFunction f);
    }
    
    public static class Tools {
    	public static void cleanIfNeeded(Minimizer minimizer) {
    		if (minimizer != null && minimizer instanceof Minimizer.NeedsCleanup) {
    			((Minimizer.NeedsCleanup)minimizer).cleanWithoutCrashing();
    		}
    	}
    }
}
