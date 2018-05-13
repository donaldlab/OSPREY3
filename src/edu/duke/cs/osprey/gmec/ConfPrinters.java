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

import java.util.Arrays;
import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

public class ConfPrinters implements ConfPrinter {
	
	public final List<ConfPrinter> printers;
	
	public ConfPrinters(ConfPrinter ... printers) {
		this.printers = Arrays.asList(printers);
	}

	@Override
	public void print(EnergiedConf conf, SimpleConfSpace confSpace, EnergyRange window) {
		for (ConfPrinter printer : printers) {
			printer.print(conf, confSpace, window);
		}
	}
}
