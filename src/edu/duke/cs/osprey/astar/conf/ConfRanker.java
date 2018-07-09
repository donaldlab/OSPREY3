/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.astar.conf;

import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;
import edu.duke.cs.osprey.tools.TimeFormatter;

import java.math.BigInteger;
import java.util.*;


public class ConfRanker {

	// TODO: use tuples pruning (e.g. pruning matrix)?

	public static interface ScorerFactory {
		AStarScorer make(EnergyMatrix emat);
	}

	public static class Builder {

		private final SimpleConfSpace confSpace;
		private final EnergyMatrix emat;

		private RCs rcs = null;
		private ScorerFactory gscorerFactory = null;
		private ScorerFactory hscorerFactory = null;
		private boolean reportProgress = false;

		public Builder(SimpleConfSpace confSpace, EnergyMatrix emat) {
			this.confSpace = confSpace;
			this.emat = emat;
		}

		public Builder setRCs(RCs val) {
			rcs = val;
			return this;
		}

		public Builder setGScorerFactory(ScorerFactory val) {
			gscorerFactory = val;
			return this;
		}

		public Builder setHScorerFactory(ScorerFactory val) {
			hscorerFactory = val;
			return this;
		}

		public Builder setReportProgress(boolean val) {
			reportProgress = val;
			return this;
		}

		public ConfRanker build() {

			if (rcs == null) {
				rcs = new RCs(confSpace);
			}

			if (gscorerFactory == null) {
				gscorerFactory = (emat) -> new PairwiseGScorer(emat);
			}

			if (hscorerFactory == null) {
				hscorerFactory = (emat) -> new TraditionalPairwiseHScorer(emat, rcs);
			}

			return new ConfRanker(confSpace, emat, rcs, gscorerFactory, hscorerFactory, reportProgress);
		}
	}

	private static int Unassigned = -1;

	private class Progress {

		public final BigInteger total;
		public BigInteger below = BigInteger.ZERO;
		public BigInteger above = BigInteger.ZERO;

		public final long startTimeNs = System.nanoTime();
		public long reportIntervalMs = 5000;

		private long lastReportTimeMs = 0;

		public Progress(BigInteger total) {
			this.total = total;
		}

		public void incrementBelow() {
			incrementBelow(BigInteger.ONE);
		}

		public void incrementAbove() {
			incrementAbove(BigInteger.ONE);
		}

		public void incrementBelow(BigInteger val) {
			below = below.add(val);
		}

		public void incrementAbove(BigInteger val) {
			above = above.add(val);
		}

		public void writeReportIfNeeded() {

			// should we write a report?
			if (reportProgress) {
				long timeMs = System.currentTimeMillis();
				if (timeMs - lastReportTimeMs >= reportIntervalMs) {
					lastReportTimeMs = timeMs;
					System.out.println(getReport());
				}
			}
		}

		public String getReport() {
			return String.format("progress: [%e,%e] of %e  (%.6f%%)   %s",
				below.doubleValue(),
				total.doubleValue() - above.doubleValue(),
				total.doubleValue(),
				below.add(above).doubleValue()/total.doubleValue()*100.0,
				TimeFormatter.format(System.nanoTime() - startTimeNs, 2)
			);
		}
	}

	private static class Node implements ConfAStarNode {

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
			Conf.index(assignments, confIndex);
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

	public final SimpleConfSpace confSpace;
	public final EnergyMatrix emat;
	public final RCs rcs;
	public final ScorerFactory gscorerFactory;
	public final ScorerFactory hscorerFactory;
	public final boolean reportProgress;

	private final AStarScorer gscorer;
	private final AStarScorer hscorer;
	private final AStarScorer negatedHScorer;
	private final ConfIndex confIndex;
	private final Node rootNode;

	private ConfRanker(SimpleConfSpace confSpace, EnergyMatrix emat, RCs rcs, ScorerFactory gscorerFactory, ScorerFactory hscorerFactory, boolean reportProgress) {

		this.confSpace = confSpace;
		this.rcs = rcs;
		this.emat = emat;
		this.gscorerFactory = gscorerFactory;
		this.hscorerFactory = hscorerFactory;
		this.reportProgress = reportProgress;

		// make the A* scorers
		gscorer = gscorerFactory.make(emat);
		hscorer = hscorerFactory.make(emat);
		negatedHScorer = hscorerFactory.make(new NegatedEnergyMatrix(confSpace, emat));

		confIndex = new ConfIndex(confSpace.positions.size());

		// make the root node
		rootNode = new Node(confSpace.positions.size());
		rootNode.index(confIndex);
		rootNode.gscore = gscorer.calc(confIndex, rcs);
		rootNode.minHScore = hscorer.calc(confIndex, rcs);
		rootNode.maxHScore = -negatedHScorer.calc(confIndex, rcs);
	}

	public BigInteger getNumConfsAtMost(double queryScore) {
		Progress progress = new Progress(rcs.getNumConformations());
		numConfsAtMost(rootNode, queryScore, progress);
		return progress.below;
	}

	private void numConfsAtMost(Node node, double queryScore, Progress progress) {

		node.index(confIndex);
		assert (confIndex.numUndefined > 0);

		// if this is the last assignment, just count the energies quickly without sub-tree bounding
		if (confIndex.numUndefined == 1) {
			countLeaves(queryScore, progress);
		} else {
			countBranches(node, queryScore, progress);
		}
	}

	private void countLeaves(double queryScore, Progress progress) {

		int pos = confIndex.undefinedPos[0];

		for (int rc : rcs.get(pos)) {

			double score = gscorer.calcDifferential(confIndex, rcs, pos, rc);

			if (score <= queryScore) {
				progress.incrementBelow();
			} else {
				progress.incrementAbove();
			}
		}

		progress.writeReportIfNeeded();
	}

	private void countBranches(Node node, double queryScore, Progress progress) {

		// find the possible assignment that maximizes the number of pruned confs
		List<Node> childNodes = new ArrayList<>();
		double bestPosScore = Double.NEGATIVE_INFINITY;
		int bestPos = -1;
		for (int i=0; i<confIndex.numUndefined; i++) {
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
				progress.incrementAbove(childNode.getNumConformations(rcs));
				iter.remove();
				continue;
			}

			if (childNode.getMaxScore() <= queryScore) {
				progress.incrementBelow(childNode.getNumConformations(rcs));
				iter.remove();
				continue;
			}

			// can't prune, keep this child node in the list
		}

		progress.writeReportIfNeeded();

		// recurse on the child nodes
		for (Node childNode: childNodes) {
			numConfsAtMost(childNode, queryScore, progress);
		}
	}
}
