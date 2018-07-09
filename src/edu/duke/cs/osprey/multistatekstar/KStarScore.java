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
