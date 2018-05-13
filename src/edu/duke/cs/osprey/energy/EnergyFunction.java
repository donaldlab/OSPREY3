/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
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

