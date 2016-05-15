/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar;

/**
 *
 * @author mhall44
 */
public interface AStarNode extends Comparable<AStarNode> {

	int[] getNodeAssignments();
	double getScore();
	void setScore(double val);
	boolean scoreNeedsRefinement();
	void setScoreNeedsRefinement(boolean val);
	int getLevel();
	boolean isFullyDefined();

	public static interface Factory<T extends AStarNode> {
		T makeRoot();
		T make(T parent, int pos, int rc);
	}
}
