package edu.duke.cs.osprey.astar.conf;

import java.math.BigInteger;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.astar.AStarProgress;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.ConfSearch;

public class ConfAStarTree implements ConfSearch {
	
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
		progress = new AStarProgress(rcs.getNumPos(), rcs.numNonTrivialPos());
	}
	
	public void stopProgress() {
		progress = null;
	}
	
	@Override
	public BigInteger getNumConformations() {
    	BigInteger num = BigInteger.valueOf(1);
    	for (int pos=0; pos<rcs.getNumPos(); pos++) {
    		num = num.multiply(BigInteger.valueOf(rcs.get(pos).length));
    	}
    	return num;
	}

	@Override
	public int[] nextConf() {
		int[] conf = new int[rcs.getNumPos()];
		nextLeafNode().getConf(conf);
		return conf;
	}
	
	public ConfAStarNode nextLeafNode() {
		
		// do we have a root node yet?
		if (rootNode == null) {
			
			// make it
			rootNode = new ConfAStarNode();
			scoreNode(rootNode);
			queue.add(rootNode);
			
			// TODO: and init progress reporting while we're at it
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
				
				// report final progress for the first leaf node, then stop reporting
				// the rest of the nodes are relatively trivial to compute
				if (progress != null) {
					progress.printProgressReport();
					progress = null;
				}
				
				return node;
			}
			
			// which pos to expand next?
			int numChildren = 0;
			confIndex.index(node);
			int nextPos = order.getNextPos(confIndex, rcs);
			for (int rc : rcs.get(nextPos)) {
				
				ConfAStarNode child = new ConfAStarNode(node, nextPos, rc);
				scoreNodeDifferential(node, child, nextPos, rc);
				
				// impossible node? skip it
				if (child.getScore() == Double.POSITIVE_INFINITY) {
					continue;
				}
				
				queue.add(child);
				numChildren++;
			}
			
            if (progress != null) {
            	progress.reportNode(node.getLevel(), node.getGScore(), node.getHScore(), queue.size(), numChildren);
            }
		}
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
