/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.energy;

import java.io.Serializable;
import java.util.List;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.AutoCleanable;

/**
 *
 * @author mhall44
 */
public interface EnergyFunction extends Serializable, AutoCleanable {
    
    public abstract double getEnergy();
    
    //we'll let each EnergyFunction object be an actual function from molecular structure to a real number (energy)
    //"partial computations" will be handled using different EnergyFunction objects

    public static interface DecomposableByDof extends EnergyFunction {
        List<EnergyFunction> decomposeByDof(Molecule m, List<DegreeOfFreedom> dofs);
    }
    
    public static interface NeedsInit extends EnergyFunction {
    	void init(Molecule m, List<DegreeOfFreedom> dofs, DoubleMatrix1D initialX);
    }
    
    public static interface NeedsCleanup extends EnergyFunction, AutoCleanable {}
    
    public static interface ExplicitChemicalChanges extends EnergyFunction {
    	int handleChemicalChanges();
    }
    
    public static class Tools {
    	
    	public static void cleanIfNeeded(EnergyFunction efunc) {
    		if (efunc != null && efunc instanceof EnergyFunction.NeedsCleanup) {
    			((EnergyFunction.NeedsCleanup)efunc).cleanWithoutCrashing();
    		}
    	}
    }

    @Override
    default void clean() {
    	Tools.cleanIfNeeded(this);
	}
}
