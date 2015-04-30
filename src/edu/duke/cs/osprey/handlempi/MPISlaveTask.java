/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.handlempi;

/**
 *
 * @author mhall44
 */
//Jobs for MPI are handled by creating a bunch of MPISlaveTask objects
//these are farmed out to various slave nodes
//they doCalculation and return whatever they're supposed to
//which can then be cast by the master to the form used by the task calling the master

public interface MPISlaveTask {
    
    Object doCalculation();
    
}
