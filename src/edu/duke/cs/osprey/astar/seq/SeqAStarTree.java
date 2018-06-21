package edu.duke.cs.osprey.astar.seq;

import edu.duke.cs.osprey.astar.ewakstar.EWAKStarLimitedSequenceTrie;
import edu.duke.cs.osprey.astar.seq.nodes.LinkedSeqAStarNode;
import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;
import edu.duke.cs.osprey.astar.seq.order.SeqAStarOrder;
import edu.duke.cs.osprey.astar.seq.scoring.SeqAStarScorer;
import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.newEwakstar.EwakstarLimitedSequenceTrie;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigInteger;
import java.util.List;


public class SeqAStarTree {

	public static class Builder {

		private final RTs rts;

		private MathTools.Optimizer optimizer = MathTools.Optimizer.Minimize;
		private SeqAStarOrder order;
		private SeqAStarScorer gscorer;
		private SeqAStarScorer hscorer;
		private int numMutable;
		private String mutableType = "max"; //defaults to max, so how anyone else was using this should be fine.
		private EwakstarLimitedSequenceTrie elst = null;

		public Builder(RTs rts) {
			this.rts = rts;
		}

		//for limiting the sequence space to a specific set of sequences in ewakstar
		public Builder setSeqTrie(EwakstarLimitedSequenceTrie obj){
			elst = obj;
			return this;
		}

		//allows for using "max", "exact", or "all" setting for how the mutable residue work.
		public Builder setMutableType(String type){
			mutableType = type;
			return this;
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

		public Builder setNumMutable(int val) {
			this.numMutable = val;
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
				mutableType,
				numMutable,
				elst
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
	public final int numMutable;
	public final String mutableType; //defaults to "max" in the above Builder
	public final EwakstarLimitedSequenceTrie elst;

	private SeqAStarNode trivialRootNode = null;
	private SeqAStarNode.Assignments assignments;

	private SeqAStarTree(RTs rts, MathTools.Optimizer optimizer, Queue<SeqAStarNode> queue, SeqAStarNode rootNode, SeqAStarOrder order, SeqAStarScorer gscorer, SeqAStarScorer hscorer, String mutableType, int numMutable, EwakstarLimitedSequenceTrie elst) {
		this.rts = rts;
		this.optimizer = optimizer;
		this.queue = queue;
		this.rootNode = rootNode;
		this.order = order;
		this.gscorer = gscorer;
		this.hscorer = hscorer;
		this.mutableType = mutableType;
		this.numMutable = numMutable;
		this.elst = elst;

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

			node.getAssignments(assignments);

			//checks if sequence is in the trie if you are using one
			boolean keepSeq = true;
			if(elst!=null && node.getLevel()!=0) {
				String seq = node.makeSequence(elst.seqSpace).toString();
				System.out.println(seq);
				if (!elst.containsSeq(seq))
					keepSeq = false;
			}



			if(keepSeq) {
				// leaf node? if mutableType = "exact", check if we have the desired number of mutations before reporting it.
				if (mutableType.equals("exact") && node.getLevel() == rts.numPos) {
					if (rts.getNumMutations(assignments) == numMutable || rts.getNumMutations(assignments) == 0) //keep wild-type sequence in your search
						return node;
					else
						continue;

				} else if (node.getLevel() == rts.numPos) {
					return node;
				}

				// which pos to expand next?
				int nextPos = order.getNextPos(assignments, rts);

				// are more mutations allowed here?
				boolean moreMutationsAllowed = rts.getNumMutations(assignments) < numMutable;

				if (moreMutationsAllowed) {

					// add all possible mutations
					for (int nextRt : rts.indicesAt(nextPos)) {
						addChild(node, nextPos, nextRt);
					}

				} else {

					// just add the wild-type at the next pos
					addChild(node, nextPos, rts.wildTypeAt(nextPos));
				}
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
