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

package edu.duke.cs.osprey.astar.seq.scoring;


import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;

public interface SeqAStarScorer {

	/**
	 * calculate a score for the given node
	 */
	double calc(SeqAStarNode.Assignments assignments);

	/**
	 * calculate a score for a given assignment, using the parent node to optimize, if possible
	 */
	default double calcDifferential(SeqAStarNode.Assignments assignments, int nextPos, int nextRt) {

		// by default update the assignments and punt to calc()
		assignments.assign(nextPos, nextRt);
		double score = calc(assignments);
		assignments.unassign(nextPos);
		return score;
	}
}
