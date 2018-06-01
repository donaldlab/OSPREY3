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

package edu.duke.cs.osprey.multistatekstar;

import java.util.Comparator;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;

/**
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */

public class ConfComparator implements Comparator<ScoredConf> {

	@Override
	public int compare(ScoredConf o1, ScoredConf o2) {
		double e1 = o1 instanceof EnergiedConf ? ((EnergiedConf)o1).getEnergy() : o1.getScore();
		double e2 = o2 instanceof EnergiedConf ? ((EnergiedConf)o2).getEnergy() : o2.getScore();
		return e1 <= e2 ? 1 : -1;
	}

}
