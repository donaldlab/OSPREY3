package edu.duke.cs.osprey.astar;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static edu.duke.cs.osprey.tools.Log.log;


public class SMAStarLab {

	public static void main(String[] args) {

		Molecule mol = PDBIO.readResource("/1CC8.ss.pdb");
		Strand strand = new Strand.Builder(mol).build();
		//List<String> resNums = Arrays.asList("A2", "A3", "A4");
		List<String> resNums = Arrays.asList("A2", "A3", "A4", "A5", "A6", "A7");
		for (String resNum : resNums) {
			strand.flexibility.get(resNum).setLibraryRotamers("VAL");
		}

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			 .setParallelism(Parallelism.makeCpu(8))
			 .build()
		) {

			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();

			RCs rcs = new RCs(confSpace);

			// enumerate the confs using A*
			ConfAStarTree astar = new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build();
			Stopwatch astarStopwatch = new Stopwatch().start();
			List<ConfSearch.ScoredConf> astarConfs = astar.nextConfs(Double.POSITIVE_INFINITY);
			astarStopwatch.stop();

			// enumerate the confs using SMA*
			SMAStar smastar = new SMAStar(emat, rcs, rcs.getNumPos() + 1);
			Stopwatch smastarStopwatch = new Stopwatch().start();
			List<ConfSearch.ScoredConf> smastarConfs = smastar.nextConfs();
			smastarStopwatch.stop();

			dumpConfs("A*", astarConfs);
			dumpConfs("SMA*", smastarConfs);

			log("  A* finished in %s", astarStopwatch.getTime(2));
			log("SMA* finished in %s", smastarStopwatch.getTime(2));

			checkConfs(astarConfs, smastarConfs);
		}
	}

	private static void dumpConfs(String name, List<ConfSearch.ScoredConf> confs) {
		log("%s confs: %d", name, confs.size());
		int size = Math.min(100, confs.size());
		for (int i=0; i<size; i++) {
			ConfSearch.ScoredConf conf = confs.get(i);
			log("\t%14.4f   %s", conf.getScore(), Conf.toString(conf.getAssignments()));
		}
	}

	private static void checkConfs(List<ConfSearch.ScoredConf> expectedConfs, List<ConfSearch.ScoredConf> observedConfs) {

		assertThat(observedConfs.size(), is(expectedConfs.size()));

		for (int i=0; i<expectedConfs.size(); i++) {
			ConfSearch.ScoredConf exp = expectedConfs.get(i);
			ConfSearch.ScoredConf obs = observedConfs.get(i);
			assertThat(obs.getAssignments(), is(exp.getAssignments()));
			assertThat(obs.getScore(), isAbsolutely(exp.getScore(), 1e-10));
		}
	}

	public static class SMAStar {

		static enum State {
			Unspawned, // never spawned yet
			Spawned, // in the child array
			Finished, // not spawned now, never need to spawn again
			Forgotten // not spawned now, might need to spawn again
		}

		static class Node implements ConfAStarNode {

			Node parent = null;
			int index = -1;
			Node[] spawnedChildren = null;
			State[] childStates = null;
			double[] forgottenScores = null;
			int depth = 0;
			double gscore = 0.0;
			double hscore = 0.0;
			double fscore = Double.NaN;
			int pos = -1;
			int rc = -1;

			public void allocateChildren(int numChildren) {
				spawnedChildren = new Node[numChildren];
				Arrays.fill(spawnedChildren, null);
				childStates = new State[numChildren];
				Arrays.fill(childStates, State.Unspawned);
				forgottenScores = new double[numChildren];
				Arrays.fill(forgottenScores, Double.NaN);
			}

			@Override
			public void index(ConfIndex index) {
				index.numDefined = depth;
				Node node = this;
				while (node.depth > 0) {
					int i = node.depth - 1;
					index.definedPos[i] = node.pos;
					index.definedRCs[i] = node.rc;
					node = node.parent;
				}
				index.node = this;
				index.sortDefined();
				index.updateUndefined();
			}

			public Node spawnChild(int pos, int rc, int index, double gscore, double hscore) {
				Node child = new Node();
				child.parent = this;
				child.index = index;
				child.depth = depth + 1;
				child.gscore = gscore;
				child.hscore = hscore;
				child.fscore = gscore + hscore;
				child.pos = pos;
				child.rc = rc;
				spawnedChildren[index] = child;

				// reset the forgotten score if needed
				if (childStates[index] == State.Forgotten) {
					forgottenScores[index] = Double.NaN;
				}

				childStates[index] = State.Spawned;
				return child;
			}

			public boolean hasSpawnableChildren() {
				for (State state : childStates) {
					if (state == State.Unspawned || state == State.Forgotten) {
						return true;
					}
				}
				return false;
			}

			public int getNextIndex(int numChildren) {

				if (spawnedChildren == null) {
					allocateChildren(numChildren);
				}

				// pick children we haven't spawned yet first
				for (int i=0; i<childStates.length; i++) {
					if (childStates[i] == State.Unspawned) {
						return i;
					}
				}

				// otherwise, pick the lowest child we've forgotten about
				int bestIndex = -1;
				double bestForgottenScore = Double.POSITIVE_INFINITY;
				for (int i=0; i<childStates.length; i++) {
					if (childStates[i] == State.Forgotten) {
						if (forgottenScores[i] < bestForgottenScore) {
							bestForgottenScore = forgottenScores[i];
							bestIndex = i;
						}
					}
				}

				if (bestIndex >= 0) {
					return bestIndex;
				}

				throw new Error("No more children to spawn");
			}

			private void removeFromParent() {

				assert (parent.spawnedChildren[index] == this);

				// TEMP
				//log("remove from parent %s", this);

				parent.spawnedChildren[index] = null;
				this.parent = null;
			}

			public String getPath() {

				assert (depth == 0 || parent != null) : String.format("%d:%d <- null", pos, rc);

				StringBuilder buf = new StringBuilder();
				Node node = this;
				if (node.depth == 0) {
					buf.append("(root)");
				}
				while (node.depth > 0) {
					if (buf.length() > 0) {
						buf.append(" <- ");
					}
					buf.append(node.pos);
					buf.append(":");
					buf.append(node.rc);
					assert (node.depth == 0 || node.parent != null) : buf.toString() + " <- null";
					node = node.parent;
				}
				return buf.toString();
			}

			@Override
			public String toString() {
				return String.format("%s   %9.4f (%9.4f)    children %s",
					getPath(),
					fscore,
					forgottenScores == null ? Double.NaN : Arrays.stream(forgottenScores)
						.filter(score -> !Double.isNaN(score))
						.min()
						.orElse(Double.NaN),
					childStates == null ? "(none)" : IntStream.range(0, childStates.length)
						.mapToObj(i -> String.format("%d:%s:%.4f:%.4f", i, childStates[i], forgottenScores[i], spawnedChildren[i] == null ? Double.NaN : spawnedChildren[i].fscore))
						.collect(Collectors.toList())
				);
			}

			@Override
			public ConfAStarNode assign(int pos, int rc) {
				throw new UnsupportedOperationException();
			}

			@Override
			public void getConf(int[] conf) {
				Node node = this;
				while (node.depth > 0) {
					conf[node.pos] = node.rc;
					node = node.parent;
				}
			}

			@Override
			public double getScore() {
				return fscore;
			}

			@Override
			public double getGScore() {
				return gscore;
			}

			@Override
			public void setGScore(double val) {
				throw new UnsupportedOperationException();
			}

			@Override
			public double getHScore() {
				return hscore;
			}

			@Override
			public void setHScore(double val) {
				throw new UnsupportedOperationException();
			}

			@Override
			public int getLevel() {
				throw new UnsupportedOperationException();
			}

			public boolean hasUnspawnedChildren() {
				for (State state : childStates) {
					if (state == State.Unspawned) {
						return true;
					}
				}
				return false;
			}

			public void backup(Queue q, int maxDepth) {

				if (hasUnspawnedChildren()) {
					return;
				}

				double oldScore = fscore;

				// update the fscore to the min of any known successor score
				fscore = Double.POSITIVE_INFINITY;
				for (Node child : spawnedChildren) {
					if (child != null) {
						fscore = Math.min(fscore, child.fscore);
					}
				}
				for (double score : forgottenScores) {
					if (!Double.isNaN(score)) {
						fscore = Math.min(fscore, score);
					}
				}

				// TEMP
				log("\tbackup   %s   %.4f -> %.4f   (%.4f)   %s", getPath(), oldScore, fscore, getMinForgottenScore(), this);
				assert (fscore >= oldScore) : "backup caused score decrease";

				// make sure scores are monotone
				for (Node child : spawnedChildren) {
					if (child != null) {
						assert (fscore <= child.fscore) : "scores must be monotone";
					}
				}
				for (double score : forgottenScores) {
					if (!Double.isNaN(score)) {
						assert (fscore <= score) : "scores must be monotone";
					}
				}

				if (fscore != oldScore) {

					// recurse if possible
					if (parent != null) {
						parent.backup(q, maxDepth);
					}
				}
			}

			public double getMinForgottenScore() {
				if (forgottenScores == null) {
					return Double.NaN;
				}
				double minScore = Double.NaN;
				for (double score : forgottenScores) {
					if (!Double.isNaN(score)) {
						if (Double.isNaN(minScore) || score < minScore) {
							minScore = score;
						}
					}
				}
				return minScore;
			}

			public boolean allChildrenSpawned() {
				for (State state : childStates) {
					if (state == State.Unspawned || state == State.Forgotten) {
						return false;
					}
				}
				return true;
			}

			public boolean allChildrenFinished() {
				for (State state : childStates) {
					if (state != State.Finished) {
						return false;
					}
				}
				return true;
			}

			public boolean hasSpawnedChildren() {
				if (childStates == null) {
					return false;
				}
				for (State state : childStates) {
					if (state == State.Spawned) {
						return true;
					}
				}
				return false;
			}
		}

		/**
		 * naive implementation for now
		 *
		 * optimize with tree-of-trees later
		 */
		static class Queue implements Iterable<Node> {

			private List<Node> nodes = new ArrayList<>();

			public void add(Node node) {

				assert (!nodes.contains(node));
				assert (Double.isFinite(node.fscore));

				nodes.add(node);
			}

			public void remove(Node node) {

				// it should be there
				assert (nodes.contains(node));

				nodes.remove(node);
			}

			public boolean contains(Node node) {
				return nodes.contains(node);
			}

			@Override
			public Iterator<Node> iterator() {
				return nodes.iterator();
			}

			public Node getLowestDeepest() {

				assert (!nodes.isEmpty());

				// first, find the lowest nodes
				double lowestScore = nodes.stream()
					.mapToDouble(node -> node.fscore)
					.min()
					.getAsDouble();
				List<Node> lowestNodes = nodes.stream()
					.filter(node -> node.fscore == lowestScore)
					.collect(Collectors.toList());
				assert (!lowestNodes.isEmpty());

				// then find the deepest node
				int deepestDepth = lowestNodes.stream()
					.mapToInt(node -> node.depth)
					.max()
					.getAsInt();
				List<Node> lowestDeepestNodes = lowestNodes.stream()
					.filter(node -> node.depth == deepestDepth)
					.collect(Collectors.toList());
				// TEMP
				//assert (lowestDeepestNodes.size() == 1) : "too many options\n\t" + String.join("\n\t", lowestDeepestNodes.stream().map(node -> node.toString()).collect(Collectors.toList()));
				return lowestDeepestNodes.get(0);
			}

			public Node removeHighestShallowestLeaf() {

				assert (!nodes.isEmpty());

				// first, find the highest nodes
				double highestScore = nodes.stream()
					.mapToDouble(node -> node.fscore)
					.max()
					.getAsDouble();
				List<Node> highestNodes = nodes.stream()
					.filter(node -> node.fscore == highestScore)
					.collect(Collectors.toList());
				assert (!highestNodes.isEmpty());

				// then, find the leaf nodes
				List<Node> leafNodes = highestNodes.stream()
					.filter(node -> !node.hasSpawnedChildren())
					.collect(Collectors.toList());
				assert (!leafNodes.isEmpty()) : "no leaves\n\t" + String.join("\n\t", highestNodes.stream().map(node -> node.toString()).collect(Collectors.toList()));

				// then, find the shallowest node
				int shallowestDepth = leafNodes.stream()
					.mapToInt(node -> node.depth)
					.min()
					.getAsInt();
				List<Node> highestShallowestNodes = highestNodes.stream()
					.filter(node -> node.depth == shallowestDepth)
					.collect(Collectors.toList());
				// TEMP
				//assert (highestShallowestNodes.size() == 1);
				Node node = highestShallowestNodes.get(0);

				remove(node);

				return node;
			}

			public int size() {
				return nodes.size();
			}

			public boolean isEmpty() {
				return nodes.isEmpty();
			}

			public String dumpScores() {
				return nodes.stream()
					.sorted(Comparator.comparing(node -> node.fscore))
					.map(node -> String.format("%.4f", node.fscore))
					.collect(Collectors.toList())
					.toString();
			}
		}

		EnergyMatrix emat;
		RCs rcs;
		long maxNumNodes;

		public SMAStar(EnergyMatrix emat, RCs rcs, long maxNumNodes) {
			this.emat = emat;
			this.rcs = rcs;
			this.maxNumNodes = maxNumNodes;
		}

		public List<ConfSearch.ScoredConf> nextConfs() {

			int numPos = rcs.getNumPos();
			List<ConfSearch.ScoredConf> confs = new ArrayList<>();
			ConfIndex confIndex = new ConfIndex(numPos);

			// get A* heuristics
			AStarScorer gscorer = new PairwiseGScorer(emat);
			AStarScorer hscorer = new TraditionalPairwiseHScorer(emat, rcs);

			// start the queue with the root node
			Queue q = new Queue();
			Node rootNode = new Node();
			rootNode.index(confIndex);
			rootNode.gscore = gscorer.calc(confIndex, rcs);
			rootNode.hscore = hscorer.calc(confIndex, rcs);
			rootNode.fscore = rootNode.gscore + rootNode.hscore;
			q.add(rootNode);

			long used = 1;

			double maxLeafScore = Double.NEGATIVE_INFINITY;
			double maxLowestDeepestScore = Double.NEGATIVE_INFINITY;

			if (maxNumNodes <= numPos) {
				throw new IllegalArgumentException(String.format("SMA* needs space for at least %d nodes for this problem (i.e., numPos + 1)", numPos + 1));
			}

			while (!q.isEmpty()) {

				// TEMP
				checkMonotonicity(q);

				Node node = q.getLowestDeepest();

				// TEMP
				log("lowest deepest: %s", node);

				// TEMP
				assert (node.fscore >= maxLowestDeepestScore) : "lowest deepest went backwards";
				maxLowestDeepestScore = node.fscore;

				// is it a leaf node?
				if (node.depth == numPos) {

					// is this the first time we've seen this leaf node?
					// NOTE: due to backing up of scores and roundoff error, the fand the g+h scores might differ ever so slightly
					// hopefully the score difference between two different confs is never this small
					double scoreEquivalance = 1e-12;
					double confScore = node.gscore;
					if (Math.abs(confScore - node.fscore) < scoreEquivalance) {

						// yup, collect the full conf
						confs.add(new ConfSearch.ScoredConf(node.makeConf(numPos), confScore));
					}

					// TEMP
					log("\n\nLEAF NODE %d   %s   %.4f >= %.4f  (%.20f)\n",
						confs.size(),
						Conf.toString(node.makeConf(numPos)),
						node.fscore,
						confScore,
						node.fscore - confScore
					);

					// TEMP
					assert (node.fscore >= maxLeafScore) : "score decreased!!";
					maxLeafScore = node.fscore;

					// and make sure we never spawn this node again
					node.fscore = Double.POSITIVE_INFINITY;
					node.parent.childStates[node.index] = State.Finished;
					node.parent.forgottenScores[node.index] = Double.POSITIVE_INFINITY;

					// TEMP
					log("\tleaf parent %s    in Q? %b", node.parent, q.contains(node.parent));
					node.parent.backup(q, numPos);

					Node n = node.parent;

					// remove the node from everything
					q.remove(node);
					node.removeFromParent();
					used--;
					log("used-- %d", used);

					// finish nodes recursively
					while (n != null) {

						if (!n.allChildrenFinished()) {
							break;
						}

						// TEMP
						log("\t\tparent finished too  %s   in Q? %b", n, q.contains(n));

						Node parent = n.parent;
						if (parent != null) {

							n.fscore = Double.POSITIVE_INFINITY;
							parent.childStates[n.index] = State.Finished;
							parent.forgottenScores[n.index] = Double.POSITIVE_INFINITY;

							n.removeFromParent();
						}

						if (q.contains(n)) {
							q.remove(n);
						}
						used--;
						log("used-- %d", used);

						n = parent;
					}

					// TEMP
					checkMonotonicity(q);

					continue;
				}

				// choose next pos
				int pos = node.depth;

				// choose the next RC
				int index = node.getNextIndex(rcs.getNum(pos));
				int rc = rcs.get(pos)[index];

				// score the child
				node.index(confIndex);
				double gscore = gscorer.calcDifferential(confIndex, rcs, pos, rc);
				double hscore = hscorer.calcDifferential(confIndex, rcs, pos, rc);
				Node child = node.spawnChild(pos, rc, index, gscore, hscore);

				// but really, don't have a lower score than the parent
				child.fscore = Math.max(child.fscore, node.fscore);

				// TEMP
				log("\tspawn child %d/%d    %s", child.index, node.childStates.length, child);
				log("\t\tnode child states %s", Arrays.toString(node.childStates));

				// TEMP
				assert (child.fscore >= node.fscore) : "scores must be monotone " + node;
				checkMonotonicity(q);

				if (!node.hasUnspawnedChildren()) {
					log("\tno unspawned children");
					node.backup(q, numPos);
				}

				// TEMP
				checkMonotonicity(q);

				if (node.allChildrenSpawned()) {
					// TEMP
					log("\tall children spawned");

					// remove it from the queue
					q.remove(node);
				}

				// TEMP
				checkMonotonicity(q);

				used++;
				// TEMP
				log("\tused++ %d", used);
				if (used > maxNumNodes) {

					// TEMP
					log("\tnodes full %d/%d    scores %s", used, maxNumNodes, q.dumpScores());

					// forget
					Node highest = q.removeHighestShallowestLeaf();

					// TEMP
					log("\tforgetting %s", highest);

					// only forget leaves!!!
					assert (!highest.hasSpawnedChildren()) : "tried to forget non-leaf " + highest.getPath();

					// tell the parent to forget us
					highest.parent.childStates[highest.index] = State.Forgotten;
					highest.parent.forgottenScores[highest.index] = highest.fscore;
					// TEMP
					log("\tparent scores: %.4f (%.4f)", highest.parent.fscore, highest.parent.getMinForgottenScore());

					// add the parent back to the queue if needed
					if (!q.contains(highest.parent)) {
						q.add(highest.parent);
						// TEMP
						log("\tadded parent: %s", highest.parent);
					}

					highest.removeFromParent();

					// TEMP
					checkMonotonicity(q);

					used--;
					// TEMP
					log("\tused-- %d", used);
				}

				// TEMP
				checkMonotonicity(q);

				// add the child to the queue
				q.add(child);

				// TEMP
				checkMonotonicity(q);
			}

			return confs;
		}

		// TEMP
		private void checkMonotonicity(Queue q) {

			for (Node node : q) {

				// check monotonicity
				if (node.spawnedChildren != null) {
					for (Node child : node.spawnedChildren) {
						if (child != null) {
							assert (node.fscore <= child.fscore) : "child monotonicity broken: " + node.toString();
						}
					}
				}

				if (node.forgottenScores != null) {
					for (double score : node.forgottenScores) {
						if (!Double.isNaN(score)) {
							assert (node.fscore <= score) : "forgotten monotonicity broken " + node.toString();
						}
					}
				}

				// check states too
				if (node.childStates != null) {
					for (int i=0; i<node.childStates.length; i++) {
						switch (node.childStates[i]) {
							case Unspawned:
								assert (node.spawnedChildren[i] == null);
								assert (Double.isNaN(node.forgottenScores[i]));
							break;
							case Spawned:
								assert (node.spawnedChildren[i] != null);
								assert (Double.isNaN(node.forgottenScores[i]));
							break;
							case Forgotten:
								assert (node.spawnedChildren[i] == null);
								assert (Double.isFinite(node.forgottenScores[i]));
							break;
							case Finished:
								assert (node.spawnedChildren[i] == null);
								assert (node.forgottenScores[i] == Double.POSITIVE_INFINITY);
							break;
						}
					}
				}
			}
		}
	}
}
