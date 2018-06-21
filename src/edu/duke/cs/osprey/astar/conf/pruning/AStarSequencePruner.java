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

package edu.duke.cs.osprey.astar.conf.pruning;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.tools.PrefixTreeSet;

import java.util.ArrayList;
import java.util.List;

public class AStarSequencePruner implements AStarPruner {

	public final SimpleConfSpace confSpace;

	private PrefixTreeSet<String> prunedSequences = new PrefixTreeSet<>();

	// NOTE: this class is NOT thread-safe!
	private List<String> sequence = null;
	private int[] assignments = null;


	public AStarSequencePruner(SimpleConfSpace confSpace) {

		this.confSpace = confSpace;

		int n = confSpace.positions.size();
		sequence = new ArrayList<>(n);
		for (int i=0; i<n; i++) {
			sequence.add(null);
		}
		assignments = new int[n];
	}

	public void add(Sequence sequence) {

		for (SimpleConfSpace.Position pos : confSpace.positions) {
			this.sequence.set(pos.index, sequence.get(pos.resNum).name);
		}

		prunedSequences.add(this.sequence);
	}

	@Override
	public boolean isPruned(ConfAStarNode node) {
		node.getConf(assignments);
		return isPruned();
	}

	@Override
	public boolean isPruned(ConfAStarNode node, int nextPos, int nextRc) {
		node.getConf(assignments);
		assignments[nextPos] = nextRc;
		return isPruned();
	}

	private boolean isPruned() {

		// convert assignments into a sequence
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			int rc = assignments[pos.index];
			if (rc == -1) {
				sequence.set(pos.index, null);
			} else {
				sequence.set(pos.index, pos.resConfs.get(rc).template.name);
			}
		}

		return prunedSequences.contains(sequence);
	}
}
