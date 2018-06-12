package edu.duke.cs.osprey.astar.seq;

import edu.duke.cs.osprey.astar.seq.nodes.LinkedSeqAStarNode;
import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;
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
		private int maxSimultaneousMutations;

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

		public Builder setMaxSimultaneousMutations(int val) {
			this.maxSimultaneousMutations = val;
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
				hscorer,
				maxSimultaneousMutations
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
	public final int maxSimultaneousMutations;

	private SeqAStarNode trivialRootNode = null;
	private SeqAStarNode.Assignments assignments;

	private SeqAStarTree(RTs rts, MathTools.Optimizer optimizer, Queue<SeqAStarNode> queue, SeqAStarNode rootNode, SeqAStarOrder order, SeqAStarScorer gscorer, SeqAStarScorer hscorer, int maxSimultaneousMutations) {
		this.rts = rts;
		this.optimizer = optimizer;
		this.queue = queue;
		this.rootNode = rootNode;
		this.order = order;
		this.gscorer = gscorer;
		this.hscorer = hscorer;
		this.maxSimultaneousMutations = maxSimultaneousMutations;

		assignments = new SeqAStarNode.Assignments(rts.numPos);
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
			for (int pos = 0; pos<rts.numPos; pos++) {
				if (rts.numTypesAt(pos) == 1) {
					node = node.assign(pos, rts.indicesAt(pos)[0]);
				}
			}
			assert (node.getLevel() == rts.getNumTrivialPos());

			// score the tail node of the chain we just created
			node.getAssignments(assignments);
			node.setGScore(gscorer.calc(assignments), optimizer);
			node.setHScore(hscorer.calc(assignments), optimizer);

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
			if (node.getLevel() == rts.numPos) {
				return node;
			}

			node.getAssignments(assignments);

			// which pos to expand next?
			int nextPos = order.getNextPos(assignments, rts);

			// are more mutations allowed here?
			boolean moreMutationsAllowed = rts.getNumMutations(assignments) < maxSimultaneousMutations;

			if (moreMutationsAllowed) {

				// add all possible mutations
				for (int nextRt : rts.indicesAt(nextPos)) {
					addChild(node, nextPos, nextRt);
				}

			} else {

				// just add the wild-type at the next pos
				int wildTypeIndex = rts.wildTypeAt(nextPos);
				if (wildTypeIndex >= 0) {
					addChild(node, nextPos, wildTypeIndex);
				}

				// wild type not allowed at next pos and we already hit the mutation limit, so just drop the sequence
			}
		}
	}

	private SeqAStarNode addChild(SeqAStarNode node, int nextPos, int nextRt) {

		// score the child node differentially against the parent node
		node.getAssignments(assignments);
		double gscore = gscorer.calcDifferential(assignments, nextPos, nextRt);
		double hscore = hscorer.calcDifferential(assignments, nextPos, nextRt);

		// immediately prune children with infinite scores
		if (Double.isInfinite(gscore) || Double.isInfinite(hscore)) {
			return null;
		}

		// make the child node and add it to the queue
		SeqAStarNode child = node.assign(nextPos, nextRt);
		child.setGScore(gscore, optimizer);
		child.setHScore(hscore, optimizer);
		queue.push(child);

		return child;
	}
}
