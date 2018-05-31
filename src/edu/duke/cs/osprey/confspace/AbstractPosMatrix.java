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

package edu.duke.cs.osprey.confspace;

import java.io.Serializable;
import java.util.Iterator;

public abstract class AbstractPosMatrix<T> implements PosMatrix<T>, Serializable {

	private int numPos;
	private int numPairs;

	protected AbstractPosMatrix(SimpleConfSpace confSpace) {
		this(confSpace.positions.size());
	}

	protected AbstractPosMatrix(int numPos) {
		this.numPos = numPos;
		this.numPairs = numPos*(numPos - 1)/2;
		allocate(numPairs);
	}

	protected abstract void allocate(int numPairs);

	@Override
	public int getNumPos() {
		return numPos;
	}

	@Override
	public int getNumPairs() {
		return numPairs;
	}

	private int getIndexNoCheck(int pos1, int pos2) {
		return pos1*(pos1 - 1)/2 + pos2;
	}

	protected int getIndex(int pos1, int pos2) {
		
		// res2 should be strictly less than res1
		if (pos2 > pos1) {
			int swap = pos1;
			pos1 = pos2;
			pos2 = swap;
		} else if (pos1 == pos2) {
			throw new Error("Can't pair residue " + pos1 + " with itself");
		}
		
		return getIndexNoCheck(pos1, pos2);
	}
	
	@Override
	public void fill(T val) {
		for (int pos1=0; pos1<getNumPos(); pos1++) {
			for (int pos2=0; pos2<pos1; pos2++) {
				set(pos1, pos2, val);
			}
		}
	}
	
	@Override
	public void fill(Iterator<T> val) {
		for (int pos1=0; pos1<getNumPos(); pos1++) {
			for (int pos2=0; pos2<pos1; pos2++) {
				set(pos1, pos2, val.next());
			}
		}
	}
}
