package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.linked.LinkedConfAStarNode;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.MPLPUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;

import java.math.BigInteger;
import java.util.*;

public class MARKStarNode implements Comparable<MARKStarNode> {

    private static AStarScorer gscorer;
    private static AStarScorer hscorer;
    private static AStarScorer negatedHScorer;
    /**
     * TODO: 1. Make MARKStarNodes spawn their own Node and MARKStarNode children.
     * TODO: 2. Make MARKStarNodes compute and update bounds correctly
     */

    private double errorUpperBound;
    private double errorLowerBound;
    private double errorBound;
    private List<MARKStarNode> children; // TODO: Pick appropriate data structure
    private Node confSearchNode = null;
    private ConfIndex confSearchIndex = null;
    private RCs confSearchRCs = null;


    private MARKStarNode(Node confNode) {
        System.out.println("New shiny node!");
        confSearchNode = confNode;
        errorUpperBound = confSearchNode.maxHScore;
        errorLowerBound = confSearchNode.minHScore;
        computeErrorBounds();
    }

    public void scoreNode(Node confNode, ConfIndex confIndex, RCs rcs) {
        confNode.gscore = gscorer.calc(confIndex, rcs);
        confNode.minHScore = hscorer.calc(confIndex, rcs);
        confNode.maxHScore = -negatedHScorer.calc(confIndex, rcs);

    }


    private double computeErrorBounds() {
        errorBound = 0;
        if(children == null || children.size() < 1)
            errorBound = errorUpperBound-errorLowerBound;
        else {
            for(MARKStarNode child: getChildren())
            {
                errorBound+= child.computeErrorBounds();
            }
        }
        System.out.println("Node error bound: "+errorBound);
        return errorBound;
    }

    public void index(ConfIndex confIndex) {
        confSearchNode.index(confIndex);
    }

    public Node getConfSearchNode() {
        return confSearchNode;
    }

    public MARKStarNode makeChild(Node child) {
        MARKStarNode newChild = new MARKStarNode(child);
        if(children == null)
            children = new ArrayList<>();
        children.add(newChild);
        return newChild;
    }


    public interface ScorerFactory {
        AStarScorer make(EnergyMatrix emat);
    }



    public static MARKStarNode makeRoot(SimpleConfSpace confSpace, EnergyMatrix energyMatrix, RCs rcs,
                                        ScorerFactory gscorerFactory, ScorerFactory hscorerFactory,
                                        boolean reportProgress) {


        AStarSearchFactory aStarSearchFactory = (emat, astarrcs) -> {
            return new ConfAStarTree.Builder(emat, astarrcs)
                    .setMPLP()
                    .build();
        };

        ConfSearch AStarTree = aStarSearchFactory.make(energyMatrix, rcs);
		// make the A* scorers
		gscorer = gscorerFactory.make(energyMatrix);
		hscorer = hscorerFactory.make(energyMatrix);
		negatedHScorer = hscorerFactory.make(new NegatedEnergyMatrix(confSpace, energyMatrix));

		ConfIndex confIndex = new ConfIndex(confSpace.positions.size());

		// make the root node
		Node rootNode = new Node(confSpace.positions.size());
		rootNode.index(confIndex);
		rootNode.gscore = gscorer.calc(confIndex, rcs);
		rootNode.minHScore = hscorer.calc(confIndex, rcs);
		rootNode.maxHScore = -negatedHScorer.calc(confIndex, rcs);
		return new MARKStarNode(rootNode);
	}


    private void dummyBoundComputationCode(Node node, double queryScore, RCs rcs) {

		// find the possible assignment that maximizes the number of pruned confs
		List<Node> childNodes = new ArrayList<>();
		double bestPosScore = Double.NEGATIVE_INFINITY;
		int bestPos = -1;
        ConfIndex confIndex = null;
        for (int i = 0; i<confIndex.numUndefined; i++) {
			int pos = confIndex.undefinedPos[i];

			int[] posRCs = rcs.get(pos);
			int numSubTreesPruned = 0;

			for (int rc : posRCs) {

				Node childNode = node.assign(pos, rc); // TODO: use object pool to re-use memory?

				// approximate the optimal sub-tree min,max scores using g+h scores
				childNode.gscore = gscorer.calcDifferential(confIndex, rcs, pos, rc);
				childNode.minHScore = hscorer.calcDifferential(confIndex, rcs, pos, rc);

				if (childNode.getMinScore() > queryScore) {
					numSubTreesPruned++;
				} else {

					childNode.maxHScore = -negatedHScorer.calcDifferential(confIndex, rcs, pos, rc);
					if (childNode.getMaxScore() <= queryScore) {
						numSubTreesPruned++;
					}
				}

				childNodes.add(childNode);
			}

			// update the best pos so far
			double posScore = (double)numSubTreesPruned/posRCs.length;
			if (posScore > bestPosScore) {
				bestPosScore = posScore;
				bestPos = pos;
			}
		}
		assert (bestPos >= 0);

		// try to prune child nodes under the best pos
		Iterator<Node> iter = childNodes.iterator();
		while (iter.hasNext()) {
			Node childNode = iter.next();

			// ignore child nodes from the non-best positions
			if (childNode.pos != bestPos) {
				iter.remove();
				continue;
			}

			if (childNode.getMinScore() > queryScore) {
				iter.remove();
				continue;
			}

			if (childNode.getMaxScore() <= queryScore) {
				iter.remove();
				continue;
			}

			// can't prune, keep this child node in the list
		}

	}

    public double getErrorLowerBound(){
        return errorLowerBound;
    }

    public double getErrorUpperBound(){
        return errorUpperBound;
    }

    public Collection<MARKStarNode> getChildren(){
        if(children == null)
            generateChildren();
        return this.children;
    }

    private void generateChildren() {
    }


    @Override
    public int compareTo(MARKStarNode other){
        return Double.compare(this.errorBound,other.errorBound);

    }

    public void expand() {
    }

    public void computeBounds() {
        if(children == null || children.size() < 1) {
            errorBound = confSearchNode.getHScore();
            return;
        }
        double errorSum = 0;
        for(MARKStarNode childNode: children) {
            errorSum += childNode.getErrorBound();
        }
        System.out.println("Error bound computed: "+errorSum);
        errorBound = errorSum;
    }

    public double getErrorBound() {
        if(children == null || children.size() < 1) {
            errorBound = confSearchNode.getHScore();
            return errorBound;
        }
        double errorSum = 0;
        for(MARKStarNode childNode: children) {
            errorSum += childNode.getErrorBound();
        }
        errorBound = errorSum;
        return errorBound;
    }


    public static class Node implements ConfAStarNode {

        private static int Unassigned = -1;
        public double gscore = Double.NaN;
        public double minHScore = Double.NaN;
        public double maxHScore = Double.NaN;
        public int[] assignments;
        public int pos = Unassigned;
        public int rc = Unassigned;

        public Node(int size) {
            assignments = new int[size];
            Arrays.fill(assignments, Unassigned);
        }

        @Override
        public Node assign(int pos, int rc) {
            Node node = new Node(assignments.length);
            node.pos = pos;
            node.rc = rc;
            System.arraycopy(assignments, 0, node.assignments, 0, assignments.length);
            node.assignments[pos] = rc;
            return node;
        }

        @Override
        public double getGScore() {
            return gscore;
        }

        @Override
        public void setGScore(double val) {
            gscore = val;
        }

        public double getMinScore() {
            return gscore + minHScore;
        }

        public double getMaxScore() {
            return gscore + maxHScore;
        }

        @Override
        public double getHScore() {
            // We want it to be as small as possible since our A* implementation finds the min
            return getMinScore()-getMaxScore();
        }

        @Override
        public void setHScore(double val) {
            throw new UnsupportedOperationException();
        }

        @Override
        public int getLevel() {
            throw new UnsupportedOperationException();
        }

        @Override
        public void getConf(int[] conf) {
            throw new UnsupportedOperationException();
        }

        @Override
        public void index(ConfIndex confIndex) {
            confIndex.numDefined = 0;
            confIndex.numUndefined = 0;
            for (int pos=0; pos<assignments.length; pos++) {
                int rc = assignments[pos];
                if (rc == -1) {
                    confIndex.undefinedPos[confIndex.numUndefined] = pos;
                    confIndex.numUndefined++;
                } else {
                    confIndex.definedPos[confIndex.numDefined] = pos;
                    confIndex.definedRCs[confIndex.numDefined] = assignments[pos];
                    confIndex.numDefined++;
                }
            }
            confIndex.node = this;
        }

        public BigInteger getNumConformations(RCs rcs) {
            BigInteger numConfs = BigInteger.ONE;
            for (int pos=0; pos<assignments.length; pos++) {
                if (assignments[pos] == Unassigned) {
                    numConfs = numConfs.multiply(BigInteger.valueOf(rcs.getNum(pos)));
                }
            }
            return numConfs;
        }
    }

    public interface AStarSearchFactory {
        ConfSearch make(EnergyMatrix energyMatrix, RCs rcs);
    }
}
