package edu.duke.cs.osprey.pruning;


import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.Progress;

import java.util.*;
import java.util.concurrent.atomic.AtomicLong;

import static edu.duke.cs.osprey.tools.Log.log;


public class TransitivePruner {

	public final SimpleConfSpace confSpace;

	public TransitivePruner(SimpleConfSpace confSpace) {
		this.confSpace = confSpace;
	}

	public void pruneSingles(PruningMatrix pmat, TaskExecutor tasks) {

		// this can take a while, so track progress
		AtomicLong numSingles = new AtomicLong(0);
		pmat.forEachUnprunedSingle((pos1, rc1) -> {
			numSingles.incrementAndGet();
			return PruningMatrix.IteratorCommand.Continue;
		});
		Progress progress = new Progress(numSingles.get());

		pmat.forEachUnprunedSingle((pos1, rc1) -> {
			tasks.submit(
				() -> hasLeaf(pmat, new RCTuple(pos1, rc1)),
				(hasLeaf) -> {
					if (!hasLeaf) {
						pmat.pruneSingle(pos1, rc1);
					}
					progress.incrementProgress();
				}
			);
			return PruningMatrix.IteratorCommand.Continue;
		});
	}

	public void prunePairs(PruningMatrix pmat, TaskExecutor tasks) {

		// this can take a while, so track progress
		AtomicLong numPairs = new AtomicLong(0);
		pmat.forEachUnprunedPair((pos1, rc1, pos2, rc2) -> {
			numPairs.incrementAndGet();
			return PruningMatrix.IteratorCommand.Continue;
		});
		Progress progress = new Progress(numPairs.get());

		pmat.forEachUnprunedPair((pos1, rc1, pos2, rc2) -> {
			tasks.submit(
				() -> hasLeaf(pmat, new RCTuple(pos1, rc1, pos2, rc2)),
				(hasLeaf) -> {
					if (!hasLeaf) {
						pmat.prunePair(pos1, rc1, pos2, rc2);
					}
					progress.incrementProgress();
				}
			);
			return PruningMatrix.IteratorCommand.Continue;
		});
	}

	public void pruneTriples(PruningMatrix pmat, TaskExecutor tasks) {

		// this can take a while, so track progress
		AtomicLong numTriples = new AtomicLong(0);
		pmat.forEachUnprunedTriple((pos1, rc1, pos2, rc2, pos3, rc3) -> {
			numTriples.incrementAndGet();
			return PruningMatrix.IteratorCommand.Continue;
		});
		Progress progress = new Progress(numTriples.get());

		pmat.forEachUnprunedTriple((pos1, rc1, pos2, rc2, pos3, rc3) -> {
			tasks.submit(
				() -> hasLeaf(pmat, new RCTuple(pos1, rc1, pos2, rc2, pos3, rc3)),
				(hasLeaf) -> {
					if (!hasLeaf) {
						pmat.pruneTriple(pos1, rc1, pos2, rc2, pos3, rc3);
					}
					progress.incrementProgress();
				}
			);
			return PruningMatrix.IteratorCommand.Continue;
		});
	}

	private boolean hasLeaf(PruningMatrix pmat, RCTuple tuple) {

		// collect all the positions: assigned first, then unassigned
		List<SimpleConfSpace.Position> positions = new ArrayList<>();
		for (int pos : tuple.pos) {
			positions.add(confSpace.positions.get(pos));
		}
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			if (!positions.contains(pos)) {
				positions.add(pos);
			}
		}

		// sort unassigned positions by number of unpruned RCs,
		// so we're likely to dead-end faster
		positions.subList(tuple.size(), positions.size())
			.sort(Comparator.comparing(pos -> pos.resConfs.stream()
				.filter(rc -> !pmat.isSinglePruned(pos.index, rc.index))
				.count()
			));

		// do DFS
		Stack<ConfIndex> stack = new Stack<>();

		// start with the DFS node representing the tuple
		ConfIndex root = new ConfIndex(positions.size());
		root.numDefined = tuple.size();
		for (int i=0; i<tuple.size(); i++) {
			root.definedPos[i] = tuple.pos.get(i);
			root.definedRCs[i] = tuple.RCs.get(i);
		}
		root.sortDefined();
		root.updateUndefined();
		stack.push(root);

		while (!stack.isEmpty()) {

			ConfIndex node = stack.pop();

			// hit a leaf node? we're done here
			if (node.numDefined == positions.size()) {
				return true;
			}

			// otherwise, expand the next pos and RCs
			SimpleConfSpace.Position pos = positions.get(node.numDefined);
			for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {

				// if this child was pruned by the pruning matrix, then skip it
				if (isPruned(pmat, node, pos.index, rc.index)) {
					continue;
				}

				// otherwise, expand it
				stack.push(node.assign(pos.index, rc.index));
			}
		}

		return false;
	}

	private boolean isPruned(PruningMatrix pmat, ConfIndex confIndex, int nextPos, int nextRc) {

		// check pairs
		for (int i=0; i<confIndex.numDefined; i++) {
			int pos = confIndex.definedPos[i];
			int rc = confIndex.definedRCs[i];
			assert (pos != nextPos || rc != nextRc);
			if (pmat.getPairwise(pos, rc, nextPos, nextRc)) {
				return true;
			}
		}

		// check triples
		if (pmat.hasHigherOrderTuples()) {

			RCTuple tuple = new RCTuple(0, 0, 0, 0, 0, 0);

			for (int i1=0; i1<confIndex.numDefined; i1++) {
				int pos1 = confIndex.definedPos[i1];
				int rc1 = confIndex.definedRCs[i1];
				assert (pos1 != nextPos || rc1 != nextRc);

				for (int i2=0; i2<i1; i2++) {
					int pos2 = confIndex.definedPos[i2];
					int rc2 = confIndex.definedRCs[i2];
					assert (pos2 != nextPos || rc2 != nextRc);

					tuple.set(pos1, rc1, pos2, rc2, nextPos, nextRc);
					tuple.sortPositions();

					if (pmat.getTuple(tuple)) {
						return true;
					}
				}
			}
		}

		return false;
	}
}
