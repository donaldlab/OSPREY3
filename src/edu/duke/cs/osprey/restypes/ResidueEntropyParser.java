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

import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.StringParsing;

public class ResidueEntropyParser {

	private ResidueEntropies resent;

	public ResidueEntropyParser(ResidueEntropies resent) {
		this.resent = resent;
	}

	public void parse(String text) {
		//It is convenient to load residue entropies into a hash map, rather than
		//into template objects, because they correspond to template names
		for (String line : FileTools.parseLines(text)) {

			// skip comments
			if (line.startsWith("%")) {
				continue;
			}

			String resType = StringParsing.getToken(line,1);
			double entropy = Double.parseDouble(StringParsing.getToken(line,2));
			resent.set(resType, entropy);
		}
	}
}
