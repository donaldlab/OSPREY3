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

package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;

public class NegatedEnergyMatrix extends ProxyEnergyMatrix {

	public NegatedEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target) {
		super(confSpace, target);
	}

	@Override
	public double getConstTerm() {
		return -super.getConstTerm();
	}

	@Override
	public void setConstTerm(double val) {
		super.setConstTerm(-val);
	}

	@Override
	public Double getOneBody(int pos, int rc) {
		return -super.getOneBody(pos, rc);
	}

	@Override
	public void setOneBody(int pos, int rc, Double val) {
		super.setOneBody(pos, rc, -val);
	}

	@Override
	public Double getPairwise(int pos1, int rc1, int pos2, int rc2) {
		return -super.getPairwise(pos1, rc1, pos2, rc2);
	}

	@Override
	public void setPairwise(int pos1, int rc1, int pos2, int rc2, Double val) {
		super.setPairwise(pos1, rc1, pos2, rc2, -val);
	}
}
