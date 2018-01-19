package edu.duke.cs.osprey.astar.conf;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.astar.AStarProgress;
import edu.duke.cs.osprey.astar.conf.linked.LinkedConfAStarFactory;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.pruning.AStarPruner;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.MPLPUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.externalMemory.EMConfAStarFactory;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.tools.ObjectPool.Checkout;

public class ConfAStarTree implements ConfSearch {

	public static class Builder {
		
		/** The energy matrix to use for pairwise residue conformation energies. */
		private EnergyMatrix emat;
		
		private RCs rcs;
		private AStarOrder order = null;
		private AStarScorer gscorer = null;
		private AStarScorer hscorer = null;
		private boolean showProgress = false;
		private ConfAStarFactory factory = new LinkedConfAStarFactory();
		private AStarPruner pruner = null;
		
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
		
		public Builder setCustom(AStarOrder order, AStarScorer gscorer, AStarScorer hscorer) {
			this.order = order;
			this.gscorer = gscorer;
			this.hscorer = hscorer;
			return this;
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
		
		/**
		 * Use external memory (eg, disk, SSD, NAS) when large A* searches
		 * cannot fit in internal memory (eg, RAM).
		 * 
		 * Use {@link ExternalMemory#setInternalLimit} to set the amount of fixed internal memory
		 * and {@link ExternalMemory#setTempDir} to set the file path for external memory.
		 */
		public Builder useExternalMemory() {
			ExternalMemory.checkInternalLimitSet();
			factory = new EMConfAStarFactory();
			return this;
		}
		
		public Builder setShowProgress(boolean val) {
			showProgress = val;
			return this;
		}

		public Builder setPruner(AStarPruner val) {
			pruner = val;
			return this;
		}
		
		public ConfAStarTree build() {
			ConfAStarTree tree = new ConfAStarTree(
				order,
				gscorer,
				hscorer,
				rcs,
				factory,
				pruner
			);
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
		 * There are two options: {@link EdgeUpdater} and {@link NodeUpdater}.
		 * In practice, the EdgeUpdater seems to work best when reference energies
		 * are used. When reference energies are not used, the NodeUpdater seems
		 * to work best.
		 */
		private MPLPUpdater updater = new NodeUpdater();
		
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
	}
	
	public final AStarOrder order;
	public final AStarScorer gscorer;
	public final AStarScorer hscorer;
	public final RCs rcs;
	public final ConfAStarFactory factory;
	public final AStarPruner pruner;
	
	private final Queue<ConfAStarNode> queue;
	private final ConfIndex confIndex;
	
	private ConfAStarNode rootNode;
	private AStarProgress progress;
	private Parallelism parallelism;
	private TaskExecutor tasks;
	private ObjectPool<ScoreContext> contexts;
	
	private ConfAStarTree(AStarOrder order, AStarScorer gscorer, AStarScorer hscorer, RCs rcs, ConfAStarFactory factory, AStarPruner pruner) {
		this.order = order;
		this.gscorer = gscorer;
		this.hscorer = hscorer;
		this.rcs = rcs;
		this.factory = factory;
		this.pruner = pruner;
		
		this.queue = factory.makeQueue(rcs);
		this.confIndex = new ConfIndex(this.rcs.getNumPos());
		
		this.rootNode = null;
		this.progress = null;
		
		this.order.setScorers(this.gscorer, this.hscorer);
		
		this.contexts = new ObjectPool<>((ingored) -> {
			ScoreContext context = new ScoreContext();
			context.index = new ConfIndex(rcs.getNumPos());
			context.gscorer = gscorer.make();
			context.hscorer = hscorer.make();
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
			
			rootNode = factory.makeRootNode(rcs.getNumPos());
			
			// pick all the single-rotamer positions now, regardless of order chosen
			// if we do them first, we basically get them for free
			// so we don't have to worry about them later in the search at all
			ConfAStarNode node = rootNode;
			for (int pos=0; pos<rcs.getNumPos(); pos++) {
				if (rcs.getNum(pos) == 1) {
					node = node.assign(pos, rcs.get(pos)[0]);
				}
			}
			assert (node.getLevel() == rcs.getNumTrivialPos());
			
			// score and add the tail node of the chain we just created
			node.index(confIndex);
			node.setGScore(gscorer.calc(confIndex, rcs));
			node.setHScore(hscorer.calc(confIndex, rcs));
			queue.push(node);
		}
		
		while (true) {
			
			// no nodes left? we're done
			if (queue.isEmpty()) {
				return null;
			}
			
			// get the next node to expand
			ConfAStarNode node = queue.poll();

			// if this node was pruned dynamically, then ignore it
			if (pruner != null && pruner.isPruned(node)) {
				continue;
			}
			
			// leaf node? report it
			if (node.getLevel() == rcs.getNumPos()) {
				
				if (progress != null) {
					progress.reportLeafNode(node.getGScore(), queue.size());
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
			List<ConfAStarNode> children = new ArrayList<>();
			for (int nextRc : rcs.get(nextPos)) {
				
				if (hasPrunedPair(confIndex, nextPos, nextRc)) {
					continue;
				}

				// if this child was pruned dynamically, then don't score it
				if (pruner != null && pruner.isPruned(node, nextPos, nextRc)) {
					continue;
				}

				tasks.submit(() -> {
					
					try (Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
						ScoreContext context = checkout.get();
						
						// score the child node differentially against the parent node
						node.index(context.index);
						ConfAStarNode child = node.assign(nextPos, nextRc);
						child.setGScore(context.gscorer.calcDifferential(context.index, rcs, nextPos, nextRc));
						child.setHScore(context.hscorer.calcDifferential(context.index, rcs, nextPos, nextRc));
						return child;
					}
					
				}, (ConfAStarNode child) -> {
					
					// collect the possible children
					if (child.getScore() < Double.POSITIVE_INFINITY) {
						children.add(child);
					}
				});
			}
			tasks.waitForFinish();
			numChildren += children.size();
			queue.pushAll(children);
			
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
		
		for (int i=0; i<confIndex.numDefined; i++) {
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
