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

package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */
public interface KStarScore {

	public enum KStarScoreType {
	    Minimized,//i.e. minimization
	    PairWiseMinimized,//pw min numerator and denominator
	    MinimizedLowerBound,//discrete numerator, pw min denominator
	    MinimizedUpperBound,//pw min numerator, discrete denominator
	    Discrete,//discrete
	    DiscreteLowerBound,//discrete numerator and denominator
	    DiscreteUpperBound;//discrete numerator and denominator
	}
	
	public enum PartitionFunctionType {
		Minimized,//i.e. minimization
		Discrete,//no min; either discrete or pw min
		UpperBound;//1+epsilon on pw min
	}
	
	public MSKStarSettings getSettings();
	public BigDecimal getScore();
	public BigDecimal getLowerBoundScore();
	public BigDecimal getUpperBoundScore();
	
	public String toString();
	public void compute(int maxNumConfs);
	public boolean constrSatisfied();
	public boolean isFullyAssigned();
	public boolean isFinal();
	public boolean isFullyProcessed();
	
}
