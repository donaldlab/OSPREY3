/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.energy;

import java.io.Serializable;

/**
 *
 * @author mhall44
 */
public interface EnergyFunction extends Serializable {
    
    public abstract double getEnergy();
    
    //we'll let each EnergyFunction object be an actual function from molecular structure to a real number (energy)
    //"partial computations" will be handled using different EnergyFunction objects


}