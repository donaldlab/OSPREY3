package edu.duke.cs.osprey.coffee.bounds;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.nodedb.NodeTree;
import edu.duke.cs.osprey.tools.BigExp;


public interface Bounder {

	/**
	 * Compute the g-score of a (possibly partial) conformation.
	 * ie, the score of the assigned positions.
	 */
	BigExp g(ConfIndex index);

	/**
	 * Comptue the h-score of a partial conformation.
	 * ie, an upper bound on the additional contribution to the g-score of any possible sub-conformation.
	 */
	BigExp h(ConfIndex index, NodeTree tree);
}
