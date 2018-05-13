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

		if (sequence.confSpace != this.confSpace) {
			throw new IllegalArgumentException("sequence conf space doesn't match pruner conf space");
		}

		for (SimpleConfSpace.Position pos : confSpace.positions) {
			this.sequence.set(pos.index, sequence.get(pos));
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
