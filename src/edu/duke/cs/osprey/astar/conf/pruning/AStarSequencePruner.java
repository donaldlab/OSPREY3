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
