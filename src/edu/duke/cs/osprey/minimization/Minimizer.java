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

package edu.duke.cs.osprey.minimization;

//Interface for minimizers.  Instantiated using an ObjectiveFunction

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.tools.AutoCleanable;

/**
 *
 * @author mhall44
 */
public interface Minimizer extends AutoCleanable {
	
    public static class Result {
    
        public DoubleMatrix1D dofValues;
        public double energy;
        
        public Result(DoubleMatrix1D dofValues, double energy) {
            this.dofValues = dofValues;
            this.energy = energy;
        }
    }

    default Result minimize() {
    	return minimizeFromCenter();
	}

	Result minimizeFromCenter();
	Result minimizeFrom(DoubleMatrix1D x);

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

	@Override
	default void clean() {
		Tools.cleanIfNeeded(this);
	}
}
