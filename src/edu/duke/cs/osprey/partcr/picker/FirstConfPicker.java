package edu.duke.cs.osprey.partcr.picker;

import java.util.List;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;

public class FirstConfPicker implements ConfPicker {

	@Override
	public ConfAStarNode pick(List<ConfAStarNode> nodes) {
		return nodes.get(0);
	}
}
