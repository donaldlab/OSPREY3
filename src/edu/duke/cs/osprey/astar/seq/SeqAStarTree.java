package edu.duke.cs.osprey.astar.seq;

import edu.duke.cs.osprey.astar.seq.nodes.LinkedSeqAStarNode;
import edu.duke.cs.osprey.astar.seq.order.SeqAStarOrder;
import edu.duke.cs.osprey.astar.seq.scoring.SeqAStarScorer;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigInteger;


public class SeqAStarTree {

	public static class Builder {

		private final RTs rts;

		private MathTools.Optimizer optimizer = MathTools.Optimizer.Minimize;
		private SeqAStarOrder order;
		private SeqAStarScorer gscorer;
		private SeqAStarScorer hscorer;

		public Builder(RTs rts) {
			this.rts = rts;
		}

		public Builder setOptimizer(MathTools.Optimizer val) {
			optimizer = val;
			return this;
		}

		public Builder setHeuristics(SeqAStarOrder order, SeqAStarScorer gscorer, SeqAStarScorer hscorer) {
			this.order = order;
			this.gscorer = gscorer;
			this.hscorer = hscorer;
			return this;
		}

		public SeqAStarTree build() {

			// don't have options for these things yet,
			// so just use the only implementations that make sense
			Queue<SeqAStarNode> queue = Queue.PriorityFactory.of(null);
			SeqAStarNode rootNode = new LinkedSeqAStarNode();

			// make sure we set all the heuristics
			if (order == null) {
				throw new IllegalArgumentException("no order heuristic set");
			}
			if (gscorer == null) {
				throw new IllegalArgumentException("no g-score heuristic set");
			}
			if (hscorer == null) {
				throw new IllegalArgumentException("no h-score heuristic set");
			}

			return new SeqAStarTree(
				rts,
				optimizer,
				queue,
				rootNode,
				order,
				gscorer,
				hscorer
			);
		}
	}

	public final RTs rts;
	public final MathTools.Optimizer optimizer;
	public final Queue<SeqAStarNode> queue;
	public final SeqAStarNode rootNode;
	public final SeqAStarOrder order;
	public final SeqAStarScorer gscorer;
	public final SeqAStarScorer hscorer;

	private SeqAStarNode trivialRootNode = null;

	private SeqAStarTree(RTs rts, MathTools.Optimizer optimizer, Queue<SeqAStarNode> queue, SeqAStarNode rootNode, SeqAStarOrder order, SeqAStarScorer gscorer, SeqAStarScorer hscorer) {
		this.rts = rts;
		this.optimizer = optimizer;
		this.queue = queue;
		this.rootNode = rootNode;
		this.order = order;
		this.gscorer = gscorer;
		this.hscorer = hscorer;
	}

	public BigInteger getNumSequences() {
		return rts.getNumSequences();
	}

	/**
	 * adds an already-enumerated node back to the tree
	 *
	 * useful for continual refinement of scores on leaf nodes
	 */
	public void add(SeqAStarNode node) {
		queue.push(node);
	}

	public SeqAStarNode nextLeafNode() {

		// do we have a trivial root node yet?
		// (ie, a root-like node, except all trivial assignments have been made)
		if (trivialRootNode == null) {

			// start at the real root node
			trivialRootNode = rootNode;

			// pick all the single-restype positions now, regardless of order heuristic
			// if we do them first, we basically get them for free
			// so we don't have to worry about them later in the search at all
			SeqAStarNode node = trivialRootNode;
			for (int pos = 0; pos<rts.numMutablePos; pos++) {
				if (rts.numTypesAt(pos) == 1) {
					node = node.assign(pos, rts.indicesAt(pos)[0]);
				}
			}
			assert (node.getLevel() == rts.getNumTrivialPos());

			// score the tail node of the chain we just created
			node.setGScore(gscorer.calc(node), optimizer);
			node.setHScore(hscorer.calc(node), optimizer);

			// and add it to the A* queue
			queue.push(node);
		}

		while (true) {

			// no nodes left? we're done
			if (queue.isEmpty()) {
				return null;
			}

			// get the next node to expand
			SeqAStarNode node = queue.poll();

			// leaf node? report it
			if (node.getLevel() == rts.numMutablePos) {
				return node;
			}

			// which pos to expand next?
			int nextPos = order.getNextPos(node, rts);
			for (int nextRt : rts.indicesAt(nextPos)) {

				// score the child node differentially against the parent node
				double gscore = gscorer.calcDifferential(node, nextPos, nextRt);
				double hscore = hscorer.calcDifferential(node, nextPos, nextRt);

				// immediately prune children with infinite scores
				if (Double.isInfinite(gscore) || Double.isInfinite(hscore)) {
					continue;
				}

				// make the child node and add it to the queue
				SeqAStarNode child = node.assign(nextPos, nextRt);
				child.setGScore(gscore, optimizer);
				child.setHScore(hscore, optimizer);
				queue.push(child);
			}
		}
	}
}
