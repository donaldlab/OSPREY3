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
