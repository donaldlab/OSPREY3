package edu.duke.cs.osprey.markstar;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.collections4.queue.CircularFifoQueue;

import edu.duke.cs.osprey.astar.AStarProgress;
import edu.duke.cs.osprey.astar.conf.ConfAStarFactory;
import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
//import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
//import edu.duke.cs.osprey.astar.conf.ConfAStarTree.Builder;
//import edu.duke.cs.osprey.astar.conf.ConfAStarTree.MPLPBuilder;
//import edu.duke.cs.osprey.astar.conf.ConfAStarTree.ScoreContext;
//import edu.duke.cs.osprey.astar.conf.linked.LinkedConfAStarFactory;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.MPLPUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ConfSearch.Splitter.Stream;
//import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.externalMemory.EMConfAStarFactory;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.tools.ObjectPool.Checkout;

public class RecursiveAStarTree implements ConfSearch {
	/**
	 * Note, I don't think I can just extend the ConfAStarTree code because it uses
	 * a builder class and is not using generics. Probably not worth modifying right
	 * now, I'll copy paste
	 */
	public static class Builder {

		/** The energy matrix to use for pairwise residue conformation energies. */
		private EnergyMatrix emat;

		private RCs rcs;
		private AStarOrder order = null;
		private AStarScorer gscorer = null;
		private AStarScorer hscorer = null;
		private AStarScorer uscorer = null;
		private boolean showProgress = false;
		private RecursiveLinkedConfAStarFactory factory = new RecursiveLinkedConfAStarFactory();

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

		public Builder setCustom(AStarOrder order, AStarScorer gscorer, AStarScorer hscorer, AStarScorer uscorer) {
			this.order = order;
			this.gscorer = gscorer;
			this.hscorer = hscorer;
			this.uscorer = uscorer;
			return this;
		}

		/**
		 * Uses the traditional estimation function to guide the tree search.
		 * {@cite Leach1998 Leach, A.R. and Lemon, A.P., 1998. Exploring the
		 * conformational space of protein side chains using dead-end elimination and
		 * the A* algorithm. Proteins Structure Function and Genetics, 33(2),
		 * pp.227-239.}
		 */
		public Builder setTraditional() {
			this.order = new DynamicHMeanAStarOrder();
			this.gscorer = new PairwiseGScorer(emat);
			this.hscorer = new TraditionalPairwiseHScorer(emat, rcs);
			this.uscorer = new UpperBoundScorer(emat);
			return this;
		}

		/**
		 * Creates an A* search using a newer estimation function based on Max Product
		 * Linear Programming (MPLP). {@cite Globerson2008 Globerson, A. and Jaakkola,
		 * T.S., 2008. Fixing max-product: Convergent message passing algorithms for MAP
		 * LP-relaxations. In Advances in neural information processing systems (pp.
		 * 553-560).}
		 * 
		 * For large designs, this A* implementation can be dramatically faster than the
		 * traditional one, and often require much less memory too.
		 */
		public Builder setMPLP() {
			setMPLP(new MPLPBuilder());
			return this;
		}

		public Builder setMPLP(MPLPBuilder builder) {
			order = new StaticScoreHMeanAStarOrder();
			gscorer = new PairwiseGScorer(emat);
			hscorer = new MPLPPairwiseHScorer(builder.updater, emat, builder.numIterations,
					builder.convergenceThreshold);
			uscorer = new UpperBoundScorer(emat);
			return this;
		}

		/**
		 * Use external memory (eg, disk, SSD, NAS) when large A* searches cannot fit in
		 * internal memory (eg, RAM).
		 * 
		 * Use {@link ExternalMemory#setInternalLimit} to set the amount of fixed
		 * internal memory and {@link ExternalMemory#setTempDir} to set the file path
		 * for external memory.
		 */
		/*
		 * public Builder useExternalMemory() { ExternalMemory.checkInternalLimitSet();
		 * factory = new EMConfAStarFactory(); return this; }
		 */

		public Builder setShowProgress(boolean val) {
			showProgress = val;
			return this;
		}

		public RecursiveAStarTree build() {
			RecursiveAStarTree tree = new RecursiveAStarTree(order, gscorer, hscorer, uscorer, rcs, factory);
			if (showProgress) {
				tree.initProgress();
			}
			return tree;
		}
	}

	public static class MPLPBuilder {

		/**
		 * The type of MPLP update step to use for each iteration.
		 * 
		 * There are two options: {@link EdgeUpdater} and {@link NodeUpdater}. In
		 * practice, the EdgeUpdater seems to work best when reference energies are
		 * used. When reference energies are not used, the NodeUpdater seems to work
		 * best.
		 */
		private MPLPUpdater updater = new NodeUpdater();

		/**
		 * The number of MPLP iterations to execute on each A* node.
		 * 
		 * This value doesn't affect the accuracy of the conformation search, only the
		 * speed.
		 * 
		 * The more iterations, the more accurate the A* estimation function will be,
		 * and fewer nodes will need to be explored to reach a leaf node. The tradeoff
		 * though is increased compute time per node explored.
		 * 
		 * Generally, it's safe to start with one iteration, then experimentally try
		 * more iterations to see if it reduces the total A* search time.
		 */
		private int numIterations = 1;

		/**
		 * If the change in energy after an iteration of the estimation function is
		 * below this threshold, MPLP will stop iterating.
		 * 
		 * This value doesn't affect the accuracy of the conformation search, only the
		 * speed.
		 * 
		 * It also has no effect if the number of iterations is 1.
		 * 
		 * For a larger number of iterations, increasing this value may reduce the time
		 * spent at each node, at the cost of exploring more total nodes. Decreasing
		 * this value may increase time spent at each node, but not necessarily reduce
		 * the total number of nodes explored.
		 * 
		 * Generally, this value won't need to be adjusted for most designs. For designs
		 * with large numbers of MPLP iterations, optimizing this value may increase
		 * performance though.
		 */
		private double convergenceThreshold = 0.0001;

		public MPLPBuilder setUpdater(MPLPUpdater val) {
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

	private static class ScoreContext {
		public ConfIndex index;
		public AStarScorer gscorer;
		public AStarScorer hscorer;
		public AStarScorer uscorer;
	}

	public final AStarOrder order;
	public final AStarScorer gscorer;
	public final AStarScorer hscorer;
	public final AStarScorer uscorer;
	public final RCs rcs;
	public final RecursiveLinkedConfAStarFactory factory;

	private final Queue<RecursiveLinkedAStarNode> leafQueue;
	private final Queue<RecursiveLinkedAStarNode> treeQueue;

	private final ConfIndex confIndex;

	private RecursiveLinkedAStarNode rootNode;
	private RecursiveLinkedAStarNode curSubTree;
	private AStarProgress progress;
	private Parallelism parallelism;
	private TaskExecutor tasks;
	private ObjectPool<ScoreContext> contexts;
	
	private boolean newSubTreeFlag;
	private double epsilon; // for now just manually set here? TODO: find some way to pass eps

	private RecursiveAStarTree(AStarOrder order, AStarScorer gscorer, AStarScorer hscorer, AStarScorer uscorer, RCs rcs,
			RecursiveLinkedConfAStarFactory factory) {
		this.order = order;
		this.gscorer = gscorer;
		this.hscorer = hscorer;
		this.uscorer = uscorer;
		this.rcs = rcs;
		this.factory = factory;

		this.leafQueue = factory.makeRecQueue(rcs);
		this.treeQueue = factory.makeRecQueue(rcs);

		this.confIndex = new ConfIndex(this.rcs.getNumPos());

		this.rootNode = null;
		this.curSubTree = null;
		this.newSubTreeFlag = true;
		this.progress = null;
		
		this.epsilon = 0.95; // TODO: pass eps

		this.order.setScorers(this.gscorer, this.hscorer);

		this.contexts = new ObjectPool<>((ingored) -> {
			ScoreContext context = new ScoreContext();
			context.index = new ConfIndex(rcs.getNumPos());
			context.gscorer = gscorer.make();
			context.hscorer = hscorer.make();
			context.uscorer = uscorer.make();
			return context;
		});

		setParallelism(null);
	}

	public void initProgress() {
		progress = new AStarProgress(rcs.getNumPos());
	}

	public AStarProgress getProgress() {
		return progress;
	}

	public void stopProgress() {
		progress = null;
	}

	public void setParallelism(Parallelism val) {

		if (val == null) {
			val = Parallelism.makeCpu(1);
		}

		parallelism = val;
		tasks = parallelism.makeTaskExecutor(1000);
		contexts.allocate(parallelism.getParallelism());
	}

	@Override
	public BigInteger getNumConformations() {
		return rcs.getNumConformations();
	}

	@Override
	public ScoredConf nextConf() {
		if (curSubTree == null) {
			// could be the first call
			curSubTree = nextSubtree(epsilon);
			// or there could be no more subtrees
			if (curSubTree == null) {
				return null;
			}
		}
		RecursiveLinkedAStarNode leafNode = nextLeafNode(curSubTree);

		if (leafNode == null) {
			// we have finished the leaves of a subtree, get a new one
			newSubTreeFlag = true;
			curSubTree = nextSubtree(epsilon);
			return nextConf();
		}

		return new ScoredConf(leafNode.makeConf(rcs.getNumPos()), leafNode.getGScore());

	}

	public RecursiveLinkedAStarNode nextSubtree(double eps) {
		// IMPT: This method does not ever pop a subtree off of the subTreeQueue.
		// It only adds to that queue, and pops off of the treeQueue
		// Also, I think we still pick next node by lower bound, so don't need to mod
		// order

		// do we have a root node yet?
		if (rootNode == null) {

			// should we have one?
			if (!rcs.hasConfs()) {
				return null;
			}
			rootNode = factory.makeRecRootNode(rcs.getNumPos());

			// pick all the single-rotamer positions now, regardless of order chosen
			// if we do them first, we basically get them for free
			// so we don't have to worry about them later in the search at all
			RecursiveLinkedAStarNode node = rootNode;
			for (int pos = 0; pos < rcs.getNumPos(); pos++) {
				if (rcs.getNum(pos) == 1) {
					node = node.assign(pos, rcs.get(pos)[0]);
				}
			}
			assert (node.getLevel() == rcs.getNumTrivialPos());

			// score and add the tail node of the chain we just created
			node.index(confIndex);
			node.setGScore(gscorer.calc(confIndex, rcs));
			node.setHScore(hscorer.calc(confIndex, rcs));
			node.setUScore(uscorer.calc(confIndex, rcs));
			treeQueue.push(node);
		}

		while (true) {
			// if no nodes left we are done?
			if (treeQueue.isEmpty()) {
				return null;
			}

			RecursiveLinkedAStarNode node = (RecursiveLinkedAStarNode) treeQueue.poll();
			// TODO : Set up recursion and rest of method.

			// if upper / lower bounds are close enough return and score tree
			// condition to return a large subtree, not sure what that is yet
			/*if (calcPhi(node) >= (1 - eps)) {
				return node;
			}*/
			// if condition does not hold and we are one from leaf return
			// tree
			if (/*calcPhi(node) <= (1 - eps) && */node.getLevel() == rcs.getNumPos() - 1) {
				// TODO: return flag to denote more work needed before estimation
				return node;
			}
			// if condition does not hold, and in middle of tree, keep looking

			// which pos to expand next? Keep same as before..., adding uscore
			int numChildren = 0;
			node.index(confIndex);
			int nextPos = order.getNextPos(confIndex, rcs);
			assert (!confIndex.isDefined(nextPos));
			assert (confIndex.isUndefined(nextPos));

			// score child nodes with tasks (possibly in parallel)
			List<RecursiveLinkedAStarNode> children = new ArrayList<>();
			for (int nextRc : rcs.get(nextPos)) {

				if (hasPrunedPair(confIndex, nextPos, nextRc)) {
					continue;
				}

				tasks.submit(() -> {

					try (Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
						ScoreContext context = checkout.get();

						// score the child node differentially against the parent node
						node.index(context.index);
						RecursiveLinkedAStarNode child = node.assign(nextPos, nextRc);
						child.setGScore(context.gscorer.calcDifferential(context.index, rcs, nextPos, nextRc));
						child.setHScore(context.hscorer.calcDifferential(context.index, rcs, nextPos, nextRc));
						// TODO: check to see if this is an error, don't know syntax of calcDifferential
						child.setUScore(context.uscorer.calc(context.index, rcs));
						return child;
					}

				}, (RecursiveLinkedAStarNode child) -> {

					// collect the possible children
					if (child.getScore() < Double.POSITIVE_INFINITY) {
						children.add(child);
					}
				});
			}
			tasks.waitForFinish();
			numChildren += children.size();
			treeQueue.pushAll(children);

			if (progress != null) {
				progress.reportInternalNode(node.getLevel(), node.getGScore(), node.getHScore(), treeQueue.size(),
						numChildren);
			}
		}

	}

	public double calcPhi(RecursiveLinkedAStarNode node) {
		/*
		 * Concerns: Do I need to use a special math function for big math? Do I need to
		 * worry about point errors with exp math? Should this go inside the node
		 * itself?
		 */
		double phi = (Math.exp(-1 * (node.getGScore() + node.getUScore()))
				- Math.exp(-1 * (node.getGScore() + node.getHScore())))
				/ (Math.exp(-1 * (node.getGScore() + node.getUScore()))
						- Math.exp(-1 * (node.getGScore() + node.getHScore())));
		return phi;
	}

	public RecursiveLinkedAStarNode nextLeafNode(RecursiveLinkedAStarNode subtree) {

		// do we have a subtree root node yet?
		if (newSubTreeFlag) {

			// should we have one?
			if (!rcs.hasConfs()) {
				return null;
			}

			RecursiveLinkedAStarNode node = subtree;
			//TODO: I think this is a problem to rescore the node, remove.
			//Although it's interesting that I saw this change dramatically...
			//node.index(confIndex);
			//node.setGScore(gscorer.calc(confIndex, rcs));
			//node.setHScore(hscorer.calc(confIndex, rcs));
			//node.setUScore(uscorer.calc(confIndex, rcs));

			leafQueue.push(node);
			newSubTreeFlag = false;
		}
		
		// check to make sure we are in the right subtree
		assert (subtree == curSubTree);

		while (true) {

			// no nodes left? we're done
			if (leafQueue.isEmpty()) {
				return null;
			}

			// get the next node to expand
			RecursiveLinkedAStarNode node = leafQueue.poll();

			// leaf node? report it
			if (node.getLevel() == rcs.getNumPos()) {

				if (progress != null) {
					progress.reportLeafNode(node.getGScore(), leafQueue.size());
				}

				return node;
			}

			// which pos to expand next?
			int numChildren = 0;
			node.index(confIndex);
			int nextPos = order.getNextPos(confIndex, rcs);
			assert (!confIndex.isDefined(nextPos));
			assert (confIndex.isUndefined(nextPos));

			// score child nodes with tasks (possibly in parallel)
			List<RecursiveLinkedAStarNode> children = new ArrayList<>();
			for (int nextRc : rcs.get(nextPos)) {

				if (hasPrunedPair(confIndex, nextPos, nextRc)) {
					continue;
				}

				tasks.submit(() -> {

					try (Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
						ScoreContext context = checkout.get();

						// score the child node differentially against the parent node
						node.index(context.index);
						RecursiveLinkedAStarNode child = node.assign(nextPos, nextRc);
						child.setGScore(context.gscorer.calcDifferential(context.index, rcs, nextPos, nextRc));
						child.setHScore(context.hscorer.calcDifferential(context.index, rcs, nextPos, nextRc));
						return child;
					}

				}, (RecursiveLinkedAStarNode child) -> {

					// collect the possible children
					if (child.getScore() < Double.POSITIVE_INFINITY) {
						children.add(child);
					}
				});
			}
			tasks.waitForFinish();
			numChildren += children.size();
			leafQueue.pushAll(children);

			if (progress != null) {
				progress.reportInternalNode(node.getLevel(), node.getGScore(), node.getHScore(), leafQueue.size(),
						numChildren);
			}
		}
	}

	public List<RecursiveLinkedAStarNode> nextLeafNodes(RecursiveLinkedAStarNode subtree, double maxEnergy) {

		if (progress != null) {
			progress.setGoalScore(maxEnergy);
		}

		List<RecursiveLinkedAStarNode> nodes = new ArrayList<>();
		while (true) {

			RecursiveLinkedAStarNode node = nextLeafNode(subtree);
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
		// I suspect this might run into problems if it tries to search across multiple
		// subtrees
		List<ScoredConf> confs = new ArrayList<>();

		if (curSubTree == null) {
			// could be the first call
			curSubTree = nextSubtree(epsilon);
			// or there could be no more subtrees
			if (curSubTree == null) {
				return null;
			}
		}
		for (RecursiveLinkedAStarNode node : nextLeafNodes(curSubTree, maxEnergy)) {
			confs.add(new ScoredConf(node.makeConf(rcs.getNumPos()), node.getGScore()));
		}
		if (leafQueue.peek() == null) {
			// if we ran out of leaf nodes in the last step
			newSubTreeFlag = true;
			curSubTree = nextSubtree(epsilon);
			//recurse, adding confs from next subtree
			confs.addAll(nextConfs(maxEnergy));
		}

		return confs;
	}

	private boolean hasPrunedPair(ConfIndex confIndex, int nextPos, int nextRc) {

		// do we even have pruned pairs?
		PruningMatrix pmat = rcs.getPruneMat();
		if (pmat == null) {
			return false;
		}

		for (int i = 0; i < confIndex.numDefined; i++) {
			int pos = confIndex.definedPos[i];
			int rc = confIndex.definedRCs[i];
			assert (pos != nextPos || rc != nextRc);
			if (pmat.getPairwise(pos, rc, nextPos, nextRc)) {
				return true;
			}
		}
		return false;
	}
}
