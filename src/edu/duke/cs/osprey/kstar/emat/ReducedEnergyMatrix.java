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

package edu.duke.cs.osprey.kstar.emat;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.KSSearchProblem;

@SuppressWarnings("serial")
public class ReducedEnergyMatrix extends EnergyMatrix {

	protected KSSearchProblem sp;

	
	public ReducedEnergyMatrix(KSSearchProblem sp, EnergyMatrix emat) {
		super(emat);
		this.sp = sp;
	}
    
    
    @Override
    public Double getOneBody(int res, int index) {
    	
    	Integer pos = sp.posNums.get(res);
    	
        return super.getOneBody(pos, index);
    }
    
    
    @Override
    public Double getPairwise(int res1, int index1, int res2, int index2) {
    	
    	Integer pos1 = sp.posNums.get(res1);
		Integer pos2 = sp.posNums.get(res2);
		
		return super.getPairwise(pos1, index1, pos2, index2);
    }
    
    
    @Override
	public HigherTupleFinder<Double> getHigherOrderTerms(int res1, int index1, int res2, int index2) {

		Integer pos1 = sp.posNums.get(res1);
		Integer pos2 = sp.posNums.get(res2);

		return super.getHigherOrderTerms(pos1, index1, pos2, index2);
	}
	
	
	@Override
	public int getNumConfAtPos(int pos) {

		Integer pos1 = sp.posNums.get(pos);

		return super.getNumConfAtPos(pos1);
	}


	@Override
	public int getNumPos() {
		return sp.posNums.size();
	}
}
