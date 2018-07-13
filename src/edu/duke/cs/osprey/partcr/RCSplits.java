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

package edu.duke.cs.osprey.partcr;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;

public class RCSplits {
	
	public static class PosInfo {
	
		private List<RCInfo> infoByRC;
	
		public PosInfo(PositionConfSpace posConfSpace) {
			infoByRC = new ArrayList<>();
			for (RC rcObj : posConfSpace.RCs) {
				infoByRC.add(new RCInfo(rcObj));
			}
		}
	}
	
	public static class RCInfo {
		
		private RC parent;
		private List<RC> children;
		
		public RCInfo(RC rcObj) {
			parent = rcObj;
			children = new ArrayList<>();
		}
		
		public RC getParent() {
			return parent;
		}
		
		public boolean isSplit() {
			return !children.isEmpty();
		}
		
		public int getNumVoxels() {
			if (isSplit()) {
				return children.size();
			} else {
				return 1;
			}
		}
		
		public List<Integer> getRCs() {
			List<Integer> rcs = new ArrayList<>();
			if (isSplit()) {
				for (RC child : children) {
					rcs.add(child.RCIndex);
				}
			} else {
				rcs.add(parent.RCIndex);
			}
			return rcs;
		}
	}
	
	private List<PosInfo> infoByPos;
	private Map<RC,RCInfo> infoByRC; // yeah, it's ok to hash on instance here
	
	public RCSplits(ConfSpace confSpace) {
	
		// init info
		infoByPos = new ArrayList<>();
		infoByRC = new HashMap<>();
		for (int pos=0; pos<confSpace.numPos; pos++) {
			PosInfo posInfo = new PosInfo(confSpace.posFlex.get(pos));
			infoByPos.add(posInfo);
			for (RCInfo rcInfo : posInfo.infoByRC) {
				infoByRC.put(rcInfo.parent, rcInfo);
			}
		}
	}
	
	public RCInfo getRCInfo(int pos, int parentRC) {
		return infoByPos.get(pos).infoByRC.get(parentRC);
	}
	
	public RCInfo getRCInfo(RC rcObj) {
		return infoByRC.get(rcObj);
	}

	public void split(RC rcObj, List<RC> splitRCs) {
		
		// update the rcinfo
		RCInfo info = infoByRC.get(rcObj);
		info.children.remove(rcObj);
		info.children.addAll(splitRCs);
		
		// update the index
		infoByRC.remove(rcObj);
		for (RC splitRCObj : splitRCs) {
			infoByRC.put(splitRCObj, info);
		}
	}

	public RCs makeRCs(int[] conf) {
		int numPos = infoByPos.size();
		List<List<Integer>> rcsByPos = new ArrayList<>(numPos);
		for (int pos=0; pos<numPos; pos++) {
			rcsByPos.add(getRCInfo(pos, conf[pos]).getRCs());
		}
		return new RCs(rcsByPos);
	}
}
