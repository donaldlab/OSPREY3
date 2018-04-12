package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.linked.LinkedConfAStarNode;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;

import java.math.BigInteger;
import java.util.*;

public class MARKStarNode implements Comparable<MARKStarNode> {

    public MARKStarNode(Node rootNode) {
    }

    public MARKStarNode(LinkedConfAStarNode childNode) {
    }

    public static interface ScorerFactory {
        AStarScorer make(EnergyMatrix emat);
    }

    private double upperBound;
    private double lowerBound;
    private final AStarScorer gscorer = null;
    private final AStarScorer hscorer= null;
    private final AStarScorer negatedHScorer= null;
    private PriorityQueue<MARKStarNode> children; // TODO: Pick appropriate data structure
    /* Missing: A variable to track our state in the conformation tree */
    private ConfTreeState confTreeSearch = null;
    private LinkedConfAStarNode confSearchNode = null;
    private ConfIndex confSearchIndex = null;
    private RCs confSearchRCs = null;

    public MARKStarNode(){
        this.upperBound = Double.NaN;
        this.lowerBound = Double.NaN;
        this.children = new PriorityQueue<MARKStarNode>();
    }


    //Function to initialize a root MARKStarNode.
    public static MARKStarNode makeRoot(SearchProblem problem) {
        return null;
    }

    public static MARKStarNode makeRoot(SimpleConfSpace confSpace, EnergyMatrix energyMatrix, RCs rcs,
                                        boolean reportProgress) {

        ScorerFactory gscorerFactory = (emat) -> new PairwiseGScorer(emat);

        ScorerFactory hscorerFactory = (emat) -> new TraditionalPairwiseHScorer(emat, rcs);


		// make the A* scorers
		AStarScorer gscorer = gscorerFactory.make(energyMatrix);
		AStarScorer hscorer = hscorerFactory.make(energyMatrix);
		AStarScorer negatedHScorer = hscorerFactory.make(new NegatedEnergyMatrix(confSpace, energyMatrix));

		ConfIndex confIndex = new ConfIndex(confSpace.positions.size());

		// make the root node
		Node rootNode = new Node(confSpace.positions.size());
		rootNode.index(confIndex);
		rootNode.gscore = gscorer.calc(confIndex, rcs);
		rootNode.minHScore = hscorer.calc(confIndex, rcs);
		rootNode.maxHScore = -negatedHScorer.calc(confIndex, rcs);
		return new MARKStarNode(rootNode, gscorer, hscoe);
	}


    private double computeUpperBound() {
		// find the possible assignment that maximizes the number of pruned confs
		List<MARKStarNode> childNodes = new ArrayList<>();
		double bestPosScore = Double.NEGATIVE_INFINITY;
		int bestPos = -1;
		ConfIndex confIndex = confSearchIndex;
		RCs rcs = confSearchRCs;
		for (int i=0; i<confIndex.numUndefined; i++) {
			int pos = confIndex.undefinedPos[i];

			int[] posRCs = rcs.get(pos);
			int numSubTreesPruned = 0;

			for (int rc : posRCs) {

				LinkedConfAStarNode childNode = confSearchNode.assign(pos, rc); // TODO: use object pool to re-use memory?

				// approximate the optimal sub-tree min,max scores using g+h scores
				double localEnergy = gscorer.calcDifferential(confIndex, rcs, pos, rc);
				double localLowerbound = -negatedHScorer.calcDifferential(confIndex, rcs, pos, rc);
				double localUpperBound = hscorer.calcDifferential(confIndex, rcs, pos, rc);



				childNodes.add(new MARKStarNode(childNode));
			}

			// update the best pos so far
			double posScore = (double)numSubTreesPruned/posRCs.length;
			if (posScore > bestPosScore) {
				bestPosScore = posScore;
				bestPos = pos;
			}
		}
		assert (bestPos >= 0);



        return 0;
    }

    private double computeLowerBound(ConfTreeState confSearchState) {
        return 0;
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


    private MARKStarNode(double upperBound, double lowerBound, PriorityQueue<MARKStarNode> children){
        this.upperBound = upperBound;
        this.lowerBound = lowerBound;
        this.children = children;
    }

    public double getLowerBound(){
        return this.lowerBound;
    }

    public double getUpperBound(){
        return this.upperBound;
    }

    public Collection<MARKStarNode> getChildren(){
        return this.children;
    }

    public void updateBounds(double upperBound, double lowerBound){
        this.upperBound = upperBound;
        this.lowerBound = lowerBound;
    }


    @Override
    public int compareTo(MARKStarNode other){
        /**
         * Compares to another MARKStarNode.
         * The node with a larger difference between upper and lower bound is less.
         *
         * @param   other the {@code MARKStarNode} to compare
         * @return  1 if the other difference is smaller, -1 if the other difference
         *          is larger, and 0 otherwise.
         */
        if ((this.upperBound - this.lowerBound) < (other.upperBound - other.lowerBound)){
            return 1;
        }else if ((this.upperBound - this.lowerBound) > (other.upperBound - other.lowerBound)){
            return -1;
        }else {
            return 0;
        }
    }

    public void expand() {
    }

    public void computeBounds() {
    }


    /* Placeholder class */
    private class ConfTreeState
    {

    }

    private static class Node implements ConfAStarNode {

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
            throw new UnsupportedOperationException();
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

}
