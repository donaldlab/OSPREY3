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
