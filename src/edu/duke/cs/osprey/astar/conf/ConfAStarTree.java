package edu.duke.cs.osprey.astar.conf;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.astar.AStarProgress;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class ConfAStarTree implements ConfSearch {

	public static class Builder {
		
		/** The energy matrix to use for pairwise residue conformation energies. */
		private EnergyMatrix emat;
		
		private RCs rcs;
		private AStarOrder order;
		private AStarScorer gscorer;
		private AStarScorer hscorer;
		
		public Builder(EnergyMatrix emat, SimpleConfSpace confSpace) {
			this(emat, new RCs(confSpace));
		}
		
		public Builder(EnergyMatrix emat, PruningMatrix pmat) {
			this(emat, new RCs(pmat));
		}
		
		public Builder(EnergyMatrix emat, RCs rcs) {
			this.emat = emat;
			this.rcs = rcs;
			
			// Jeff: MPLP is dramatically faster for large A* searches
			// and for small searches, who cares how fast A* is,
			// so I think it makes a good default for all cases
			setMPLP();
		}
		
		/**
		 * Uses the traditional estimation function to guide the tree search.
		 * {@cite Leach1998 Leach, A.R. and Lemon, A.P., 1998. Exploring the conformational
		 * space of protein side chains using dead-end elimination and the A* algorithm.
		 * Proteins Structure Function and Genetics, 33(2), pp.227-239.}
		 */
		public Builder setTraditional() {
			this.order = new DynamicHMeanAStarOrder();
			this.gscorer = new PairwiseGScorer(emat);
			this.hscorer = new TraditionalPairwiseHScorer(emat, rcs);
			return this;
		}
		
		/**
		 * Creates an A* search using a newer estimation function based on Max Product Linear
		 * Programming (MPLP).
		 * {@cite Globerson2008 Globerson, A. and Jaakkola, T.S., 2008. Fixing max-product: Convergent message passing
		 * algorithms for MAP LP-relaxations. In Advances in neural information processing systems (pp. 553-560).}
		 * 
		 * For large designs, this A* implementation can be dramatically faster than the traditional
		 * one, and often require much less memory too.
		 */
		public Builder setMPLP() {
			setMPLP(new MPLPBuilder());
			return this;
		}
		
		public Builder setMPLP(MPLPBuilder builder) {
			order = new StaticScoreHMeanAStarOrder();
			gscorer = new PairwiseGScorer(emat);
			hscorer = new MPLPPairwiseHScorer(
				builder.updater,
				emat,
				builder.numIterations,
				builder.convergenceThreshold
			);
			return this;
		}
		
		public ConfSearch build() {
			return new ConfAStarTree(
				order,
				gscorer,
				hscorer,
				rcs
			);
		}
	}
	
	public static class MPLPBuilder {
		
		private NodeUpdater updater = new NodeUpdater();
		
		/**
		 * The number of MPLP iterations to execute on each A* node.
		 * 
		 * This value doesn't affect the accuracy of the conformation search, only the speed.
		 * 
		 * The more iterations, the more accurate the A* estimation function will be,
		 * and fewer nodes will need to be explored to reach a leaf node. The tradeoff though is
		 * increased compute time per node explored.
		 * 
		 * Generally, it's safe to start with one iteration, then experimentally try more
		 * iterations to see if it reduces the total A* search time.
		 */
		private int numIterations = 1;
		
		/**
		 * If the change in energy after an iteration of the estimation function is below this
		 * threshold, MPLP will stop iterating.
		 * 
		 * This value doesn't affect the accuracy of the conformation search, only the speed.
		 * 
		 * It also has no effect if the number of iterations is 1.
		 * 
		 * For a larger number of iterations, increasing this value may reduce the time spent
		 * at each node, at the cost of exploring more total nodes. Decreasing this value may
		 * increase time spent at each node, but not necessarily reduce the total number of
		 * nodes explored.
		 * 
		 * Generally, this value won't need to be adjusted for most designs. For designs with
		 * large numbers of MPLP iterations, optimizing this value may increase performance though.
		 */
		private double convergenceThreshold = 0.0001;
		
		public MPLPBuilder setUpdater(NodeUpdater val) {
			updater = val;
			return this;
		}
		
		public MPLPBuilder setNumIterations(int val) {
			numIterations = val;
			return this;
		}
		
		public MPLPBuilder setConvergenceThreshold(double val) {
			convergenceThreshold = val;
			return this;
		}
	}

	public static MPLPBuilder MPLPBuilder() {
		return new MPLPBuilder();
	}
		

	private AStarOrder order;
	private AStarScorer gscorer;
	private AStarScorer hscorer;
	private PriorityQueue<ConfAStarNode> queue;
	private RCs rcs;
	private ConfAStarNode rootNode;
	private ConfIndex confIndex;
	private AStarProgress progress;
	
	public ConfAStarTree(AStarOrder order, AStarScorer gscorer, AStarScorer hscorer, RCs rcs) {
		this.order = order;
		this.gscorer = gscorer;
		this.hscorer = hscorer;
		this.queue = new PriorityQueue<>();
		this.rcs = rcs;
		this.rootNode = null;
		this.confIndex = new ConfIndex(this.rcs.getNumPos());
		this.progress = null;
		
		this.order.setScorers(this.gscorer, this.hscorer);
	}
	
	public void initProgress() {
		progress = new AStarProgress(rcs.getNumPos());
	}
	
	public void stopProgress() {
		progress = null;
	}
	
	@Override
	public BigInteger getNumConformations() {
		
		if (rcs.hasConfs()) {
			
			BigInteger num = BigInteger.valueOf(1);
			for (int pos=0; pos<rcs.getNumPos(); pos++) {
				num = num.multiply(BigInteger.valueOf(rcs.get(pos).length));
			}
			return num;
			
		} else {
			
			return BigInteger.ZERO;
		}
	}

	@Override
	public ScoredConf nextConf() {
		ConfAStarNode leafNode = nextLeafNode();
		if (leafNode == null) {
			return null;
		}
		return new ScoredConf(
			leafNode.makeConf(rcs.getNumPos()),
			leafNode.getGScore()
		);
	}
	
	public ConfAStarNode nextLeafNode() {
		
		// do we have a root node yet?
		if (rootNode == null) {
			
			// should we have one?
			if (!rcs.hasConfs()) {
				return null;
			}
			
			rootNode = new ConfAStarNode();
			
			// pick all the single-rotamer positions now, regardless of order chosen
			// if we do them first, we basically get them for free
			// so we don't have to worry about them later in the search at all
			ConfAStarNode node = rootNode;
			for (int pos=0; pos<rcs.getNumPos(); pos++) {
				if (rcs.getNum(pos) == 1) {
					node = new ConfAStarNode(node, pos, rcs.get(pos)[0]);
				}
			}
			assert (node.getLevel() == rcs.getNumTrivialPos());
			
			// score and add the tail node of the chain we just created
			scoreNode(node);
			queue.add(node);
		}
		
		while (true) {
			
			// no nodes left? we're done
			if (queue.isEmpty()) {
				return null;
			}
			
			// get the next node to expand
			ConfAStarNode node = queue.poll();
			
			// leaf node? report it
			if (node.getLevel() == rcs.getNumPos()) {
				
				if (progress != null) {
					progress.reportLeafNode(node.getGScore());
				}
			
				return node;
			}
			
			// which pos to expand next?
			int numChildren = 0;
			confIndex.index(node);
			int nextPos = order.getNextPos(confIndex, rcs);
			assert (!confIndex.isDefined(nextPos));
			assert (confIndex.isUndefined(nextPos));
			
			for (int nextRc : rcs.get(nextPos)) {
				
				if (hasPrunedPair(confIndex, nextPos, nextRc)) {
					continue;
				}
				
				ConfAStarNode child = new ConfAStarNode(node, nextPos, nextRc);
				scoreNodeDifferential(node, child, nextPos, nextRc);
				
				// impossible node? skip it
				if (child.getScore() == Double.POSITIVE_INFINITY) {
					continue;
				}
				
				queue.add(child);
				numChildren++;
			}
			
            if (progress != null) {
            	progress.reportInternalNode(node.getLevel(), node.getGScore(), node.getHScore(), queue.size(), numChildren);
            }
		}
	}
	
	public List<ConfAStarNode> nextLeafNodes(double maxEnergy) {
		
		if (progress != null) {
			progress.setGoalScore(maxEnergy);
		}
		
		List<ConfAStarNode> nodes = new ArrayList<>();
		while (true) {
			
			ConfAStarNode node = nextLeafNode();
			if (node == null) {
				break;
			}
			
			nodes.add(node);
			
			if (node.getGScore() >= maxEnergy) {
				break;
			}
		}
		return nodes;
	}
	
	@Override
	public List<ScoredConf> nextConfs(double maxEnergy) {
		List<ScoredConf> confs = new ArrayList<>();
		for (ConfAStarNode node : nextLeafNodes(maxEnergy)) {
			confs.add(new ScoredConf(
				node.makeConf(rcs.getNumPos()),
				node.getGScore()
			));
		}
		return confs;
	}
	
	private boolean hasPrunedPair(ConfIndex confIndex, int nextPos, int nextRc) {
		
		// do we even have pruned pairs?
		PruningMatrix pmat = rcs.getPruneMat();
		if (pmat == null) {
			return false;
		}
		
		for (int i=0; i<confIndex.getNumDefined(); i++) {
			int pos = confIndex.getDefinedPos()[i];
			int rc = confIndex.getDefinedRCs()[i];
			assert (pos != nextPos || rc != nextRc);
			if (pmat.getPairwise(pos, rc, nextPos, nextRc)) {
				return true;
			}
		}
		return false;
	}

	private void scoreNode(ConfAStarNode node) {
		confIndex.index(node);
		node.setGScore(gscorer.calc(confIndex, rcs));
		node.setHScore(hscorer.calc(confIndex, rcs));
	}
	
	private void scoreNodeDifferential(ConfAStarNode parent, ConfAStarNode child, int nextPos, int nextRc) {
		confIndex.index(parent);
		child.setGScore(gscorer.calcDifferential(confIndex, rcs, nextPos, nextRc));
		child.setHScore(hscorer.calcDifferential(confIndex, rcs, nextPos, nextRc));
	}
}
