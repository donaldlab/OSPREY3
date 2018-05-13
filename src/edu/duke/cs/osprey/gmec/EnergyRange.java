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

package edu.duke.cs.osprey.gmec;

public class EnergyRange {
	
	private double min;
	private double size;
	
	public EnergyRange(double energy, double size) {
		min = energy;
		this.size = size;
	}
	
	public boolean updateMin(double energy) {
		if (energy >= min) {
			return false;
		}
		min = energy;
		return true;
	}
	
	public double getMin() {
		return min;
	}
	
	public double getMax() {
		return min + size;
	}

	public boolean contains(double energy) {
		return energy >= min && energy <= getMax();
	}
	
	public boolean containsOrBelow(double energy) {
		return energy <= getMax();
	}
}
