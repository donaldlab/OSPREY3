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

package edu.duke.cs.osprey.kstar.pruning;

import edu.duke.cs.osprey.astar.comets.UpdatedPruningMatrix;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.kstar.KSSearchProblem;

@SuppressWarnings("serial")
public class UnprunedPruningMatrix extends ReducedPruningMatrix {

	// one body, pairwise, and higher order will be called by an object that sees
	// only the reduced number of positions. when accessing RCs in the pruneMat, 
	// it is necessary to convert from reduced position to the absolute position.
	// however, reducedAllowedAAs in sp has also been cut down, so we access that using 
	// the position specified by the caller

	public UnprunedPruningMatrix(KSSearchProblem sp, UpdatedPruningMatrix upm, double pruningInterval) {
		super(sp, upm);
		this.setPruningInterval(pruningInterval);
	}


	@Override
	public Boolean getOneBody(int res, int index) {

		Integer pos = sp.posNums.get(res);
		String rcAAType = sp.confSpace.posFlex.get(pos).RCs.get(index).AAType;

		// if in specified aa list, invert
		if(sp.reducedAllowedAAs.get(res).contains(rcAAType))
			return false;

		return true;
	}


	@Override
	public Boolean getPairwise(int res1, int index1, int res2, int index2) {

		Integer pos1 = sp.posNums.get(res1);
		Integer pos2 = sp.posNums.get(res2);

		String rcAAType1 = sp.confSpace.posFlex.get(pos1).RCs.get(index1).AAType;
		String rcAAType2 = sp.confSpace.posFlex.get(pos2).RCs.get(index2).AAType;

		// if in specified aa list, invert
		if(sp.reducedAllowedAAs.get(res1).contains(rcAAType1) && sp.reducedAllowedAAs.get(res2).contains(rcAAType2))
			return false;

		return true;
	}


	@Override
	public HigherTupleFinder<Boolean> getHigherOrderTerms(int res1, int index1, int res2, int index2) {

		Integer pos1 = sp.posNums.get(res1);
		Integer pos2 = sp.posNums.get(res2);

		String rcAAType1 = sp.confSpace.posFlex.get(pos1).RCs.get(index1).AAType;
		String rcAAType2 = sp.confSpace.posFlex.get(pos2).RCs.get(index2).AAType;

		// allowedAAs has been reduced, so use resX positions
		if(sp.reducedAllowedAAs.get(res1).contains(rcAAType1) && sp.reducedAllowedAAs.get(res2).contains(rcAAType2))
			return upm.getHigherOrderTerms(pos1, index1, pos2, index2);

		return null;
	}
}
