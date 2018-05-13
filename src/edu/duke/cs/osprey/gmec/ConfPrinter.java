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
