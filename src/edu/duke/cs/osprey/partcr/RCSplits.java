package edu.duke.cs.osprey.partcr;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
		
		public boolean isSplit() {
			return !children.isEmpty();
		}
		
		public boolean isParent(int rc) {
			
			if (isSplit()) {
				throw new IllegalStateException("parent rc index is meaningless after splits");
			}
			
			return parent.RCIndex == rc;
		}
		
		public boolean isChild(int rc) {
			
			if (!isSplit()) {
				throw new IllegalStateException("no child RCs since nothing was split");
			}
			
			// TODO: could switch to constant-time lookup if we need the speed
			for (RC rcObj : children) {
				if (rcObj.RCIndex == rc) {
					return true;
				}
			}
			return false;
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
}
