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

package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

public interface ConfPrinter {
	
	void print(ConfSearch.EnergiedConf conf, SimpleConfSpace confSpace, EnergyRange range);
	
	default void print(ConfSearch.EnergiedConf conf, SimpleConfSpace confSpace) {
		print(conf, confSpace, null);
	}

	default void print(ConfSearch.ScoredConf conf, SimpleConfSpace confSpace) {
		print(new ConfSearch.EnergiedConf(conf, Double.NaN), confSpace, null);
	}
	
	default void print(ConfSearch.EnergiedConf conf) {
		print(conf, null, null);
	}

	default void print(ConfSearch.ScoredConf conf) {
		print(new ConfSearch.EnergiedConf(conf, Double.NaN), null, null);
	}

	default void cleanup() {
		// nothing to do
	}
	
	public static class Nop implements ConfPrinter {

		@Override
		public void print(ConfSearch.EnergiedConf conf, SimpleConfSpace confSpace, EnergyRange range) {
			// do nothing, ie, a no-op
		}
	}
}
