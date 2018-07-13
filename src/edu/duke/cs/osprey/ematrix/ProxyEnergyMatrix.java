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

package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

public class ProxyEnergyMatrix extends EnergyMatrix {

	public final EnergyMatrix target;

	public ProxyEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target) {
		super(confSpace);
		this.target = target;
	}

	@Override
	protected void allocate(int numOneBody, int numPairwise) {
		// don't allocate anything
		// we're just going to proxy all methods to the target emat
	}

	@Override
	public Double getOneBody(int pos, int rc) {
		return target.getOneBody(pos, rc);
	}

	@Override
	public void setOneBody(int pos, int rc, Double val) {
		target.setOneBody(pos, rc, val);
	}

	@Override
	public Double getPairwise(int pos1, int rc1, int pos2, int rc2) {
		return target.getPairwise(pos1, rc1, pos2, rc2);
	}

	@Override
	public void setPairwise(int pos1, int rc1, int pos2, int rc2, Double val) {
		target.setPairwise(pos1, rc1, pos2, rc2, val);
	}

	@Override
	public void setHigherOrder(RCTuple tup, Double val) {
		super.setHigherOrder(tup, val);
	}

	// TODO: need to proxy anything else?
}
