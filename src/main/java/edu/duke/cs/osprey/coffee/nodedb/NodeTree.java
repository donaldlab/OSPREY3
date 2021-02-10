package edu.duke.cs.osprey.coffee.nodedb;

import edu.duke.cs.osprey.astar.conf.RCs;


public class NodeTree {

	public final RCs rcs;
	public final Integer maxSimultaneousMutations;

	public NodeTree(RCs rcs, Integer maxSimultaneousMutations) {
		this.rcs = rcs;
		this.maxSimultaneousMutations = maxSimultaneousMutations;
	}

	/**
	 * Make a node tree with no restriction on the number of simultaneous mutations.
	 */
	public NodeTree(RCs rcs) {
		this(rcs, null);
	}
}
