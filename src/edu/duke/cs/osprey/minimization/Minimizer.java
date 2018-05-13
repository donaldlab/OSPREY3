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

