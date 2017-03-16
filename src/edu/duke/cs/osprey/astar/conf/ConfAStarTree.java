package edu.duke.cs.osprey.astar.conf;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.astar.AStarProgress;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class ConfAStarTree implements ConfSearch {
	
	public final AStarOrder order;
	public final AStarScorer gscorer;
	public final AStarScorer hscorer;
	public final RCs rcs;
	
	private PriorityQueue<ConfAStarNode> queue;
	private ConfAStarNode rootNode;
	private ConfIndex confIndex;
	private AStarProgress progress;
	
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
