package edu.duke.cs.osprey.partcr.picker;

import java.util.List;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;

public interface ConfPicker {
	
	ConfAStarNode pick(List<ConfAStarNode> nodes);
}
