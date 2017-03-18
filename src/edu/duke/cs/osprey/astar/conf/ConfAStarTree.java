package edu.duke.cs.osprey.astar.conf;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.astar.AStarProgress;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectPool;

public class ConfAStarTree implements ConfSearch {
	
	private static class ScoreContext {
		public ConfIndex index;
		public AStarScorer gscorer;
		public AStarScorer hscorer;
	}
	
	public final AStarOrder order;
	public final AStarScorer gscorer;
	public final AStarScorer hscorer;
	public final RCs rcs;
	
	private PriorityQueue<ConfAStarNode> queue;
	private ConfAStarNode rootNode;
	private ConfIndex confIndex;
	private AStarProgress progress;
	private Parallelism parallelism;
	private TaskExecutor tasks;
	private ObjectPool<ScoreContext> contexts;
	
	public ConfAStarTree(AStarOrder order, AStarScorer gscorer, AStarScorer hscorer, RCs rcs) {
		this.order = order;
		this.gscorer = gscorer;
		this.hscorer = hscorer;
		this.rcs = rcs;
		this.queue = new PriorityQueue<>();
		this.rootNode = null;
		this.confIndex = new ConfIndex(this.rcs.getNumPos());
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
			confIndex.index(node);
			node.setGScore(gscorer.calc(confIndex, rcs));
			node.setHScore(hscorer.calc(confIndex, rcs));
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
					progress.reportLeafNode(node.getGScore(), queue.size());
				}
			
				return node;
			}
			
			// which pos to expand next?
			int numChildren = 0;
			confIndex.index(node);
			int nextPos = order.getNextPos(confIndex, rcs);
			assert (!confIndex.isDefined(nextPos));
			assert (confIndex.isUndefined(nextPos));
			
			// score child nodes with tasks (possibly in parallel)
			List<ConfAStarNode> children = new ArrayList<>();
			for (int nextRc : rcs.get(nextPos)) {
				
				if (hasPrunedPair(confIndex, nextPos, nextRc)) {
					continue;
				}
				
				ConfAStarNode child = new ConfAStarNode(node, nextPos, nextRc);
				
				tasks.submit(() -> {
					
					try (ObjectPool<ScoreContext>.Checkout checkout = contexts.autoCheckout()) {
						ScoreContext context = checkout.get();
						
						// score the child node differentially against the parent node
						context.index.index(node);
						child.setGScore(context.gscorer.calcDifferential(context.index, rcs, nextPos, nextRc));
						child.setHScore(context.hscorer.calcDifferential(context.index, rcs, nextPos, nextRc));
					}
					
				}, (Runnable task) -> {
					
					// collect the possible children
					if (child.getScore() < Double.POSITIVE_INFINITY) {
						children.add(child);
					}
				});
			}
			tasks.waitForFinish();
			numChildren += children.size();
			queue.addAll(children);
			
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
}
