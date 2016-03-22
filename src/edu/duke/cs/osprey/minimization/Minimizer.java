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
    
    DoubleMatrix1D minimize();//return argmin for the ObjectiveFunction
    //(the minimum can then be obtained by calling the objective function on those values)
}
