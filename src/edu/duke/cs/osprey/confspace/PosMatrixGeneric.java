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
