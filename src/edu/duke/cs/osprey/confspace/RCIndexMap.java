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

package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.List;

public class RCIndexMap {
	
	// this is a map between old rc indices and new rc indices after splits
	// old rcs can map to multiple new rcs
	// new rcs can map to only one old rc
	
	private List<List<Integer>> oldToNew;
	private List<Integer> newToOld;
	
	public RCIndexMap(int numOldRCs) {
		
		// start with a simple bijection between rcs
		oldToNew = new ArrayList<>(numOldRCs);
		newToOld = new ArrayList<>(numOldRCs);
		for (int i=0; i<numOldRCs; i++) {
			oldToNew.add(makeList(i));
			newToOld.add(i);
		}
	}
	
	private List<Integer> makeList(int val) {
		List<Integer> list = new ArrayList<>();
		list.add(val);
		return list;
	}

	public void remove(int oldIndex) {
		
		// break link between old/new rcs at this index
		oldToNew.get(oldIndex).clear();
		newToOld.set(oldIndex, null);
	}

	public void add(int oldIndex, int newIndex) {
		
		// add old->new link
		oldToNew.get(oldIndex).add(newIndex);
		
		// add new->old link
		if (newIndex < newToOld.size()) {
			newToOld.set(newIndex, null);
		} else if (newIndex == newToOld.size()) {
			newToOld.add(null);
		} else {
			throw new IllegalArgumentException(String.format("can't handle new index %d yes, add %d first",
				newIndex, newToOld.size()
			));
		}
	}
	
	public Integer newToOld(int rc) {
		return newToOld.get(rc);
	}

	public List<Integer> oldToNew(int rc) {
		return oldToNew.get(rc);
	}
}
