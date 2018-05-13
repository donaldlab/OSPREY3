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

package edu.duke.cs.osprey.restypes;

import java.io.Serializable;
import java.util.HashMap;

public class ResidueEntropies implements Serializable {

	private HashMap<String,Double> values = new HashMap<>();

	public void set(String resType, double val) {
		values.put(resType, val);
	}

	public double get(String resType) {

		// return defined entropy, if any
		Double val = values.get(normalize(resType));
		if (val != null) {
			return val;
		}

		// or zero entropy, by default
		return 0;
	}

	public void clear() {
		values.clear();
	}

	public void setAll(ResidueEntropies other) {
		values.putAll(other.values);
	}

	private String normalize(String resType) {
		return resType.toUpperCase();
	}
}
