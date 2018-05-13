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
