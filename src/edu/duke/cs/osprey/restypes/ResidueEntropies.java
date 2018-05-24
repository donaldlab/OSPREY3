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
