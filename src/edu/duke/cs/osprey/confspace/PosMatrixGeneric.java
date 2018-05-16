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
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
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

public class PosMatrixGeneric<T> extends AbstractPosMatrix<T> {

	private T[] vals;

	public PosMatrixGeneric(SimpleConfSpace confSpace) {
		super(confSpace);
	}

	public PosMatrixGeneric(int numPos) {
		super(numPos);
	}

	@Override
	@SuppressWarnings("unchecked")
	protected void allocate(int numPairs) {
		vals = (T[])new Object[numPairs];
	}

	@Override
	public T get(int pos1, int pos2) {
		return vals[getIndex(pos1, pos2)];
	}

	@Override
	public void set(int pos1, int pos2, T val) {
		vals[getIndex(pos1, pos2)] = val;
	}
}




