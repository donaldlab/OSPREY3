package edu.duke.cs.osprey.markstar.visualizer;

import edu.duke.cs.osprey.tools.JvmMem;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.*;
import java.util.function.Consumer;

import static edu.duke.cs.osprey.tools.Log.log;


/**
 * Sorts the chosen levels of a conf tree in the desired order
 * and pushes those levels to the top of the tree.
 * Then outputs the chosen levels of the tree (ie, the top) to the output file.
 */
public class KStarTreeSplitter {

	public static void main(String[] args) {

		// read the args
		File inFile;
		File outFile;
		List<String> resNums;
		try {
			inFile = new File(args[0]);
			outFile = new File(args[1]);
			resNums = Arrays.asList(args[2].split(","));
		} catch (Throwable t) {
			log("Invalid arguments, caused an exception:");
			t.printStackTrace(System.out);
			log("expected arguments: inFile outFile resNums");
			log("\twhere resNums is a comma-separated list of residue numbers");
			return;
		}

		// read the tree
		Stopwatch readsw = new Stopwatch().start();
		log("reading tree from %s ...", inFile.getAbsolutePath());
		KStarTreeNode root = KStarTreeNode.parseTree(inFile, false, null);
		log("\tdone in %s, used %s", readsw.stop().getTime(2), JvmMem.getOldPool());

		log("target residue order:  %s", resNums);

		// get the current residue order
		List<String> currentResNums = getResNums(root);

		// BUBBLE SORT, BABY!!
		log("current residue order: %s", currentResNums);
		log("sorting ...");
		Stopwatch sortsw = new Stopwatch().start();
		for (int i=0; i<resNums.size(); i++) {
			for (int j=currentResNums.indexOf(resNums.get(i)); j>i; j--) {

				// swap the res nums
				int jm1 = j - 1;
				String temp = currentResNums.get(j);
				currentResNums.set(j, currentResNums.get(jm1));
				currentResNums.set(jm1, temp);

				// swap the tree layers
				pushUpTreeLayer(root, j);

				log("current residue order: %s, used %s", getResNums(root), JvmMem.getOldPool());
			}
		}
		log("sort finished in %s", sortsw.stop().getTime(2));

		Stopwatch writesw = new Stopwatch().start();
		log("writing tree to %s ...", outFile.getAbsolutePath());
		try (FileWriter writer = new FileWriter(outFile)) {
			root.printTreeLikeMARKStar(writer);
		} catch (IOException ex) {
			throw new Error(ex);
		}
		log("\tfinished in %s", writesw.stop().getTime(2));
	}

	private static List<String> getResNums(KStarTreeNode node) {

		LinkedHashSet<String> assignedResNums = new LinkedHashSet<>();
		while (true) {

			// which new residue was assigned here?
			Integer i = node.getAssignmentIndex();
			if (i != null) {
				assignedResNums.add(node.getAssignments()[i].split(":")[0]);
			}

			// go to an arbitrary child node
			if (node.children.isEmpty()) {
				break;
			}
			node = node.children.get(0);
		}

		return new ArrayList<>(assignedResNums);
	}

	private static void pushUpTreeLayer(KStarTreeNode root, int i) {

		// convert the res order index to a depth
		int depth = i + 1;
		if (depth > root.getConfAssignments().length) {
			throw new IllegalArgumentException("depth too deep");
		}
		if (depth < 2) {
			throw new IllegalArgumentException("depth too shallow");
		}

		MathContext mathContext = new MathContext(128, RoundingMode.HALF_UP);

		class Edge {

			final KStarTreeNode nodeim1;
			final int posim1;
			final KStarTreeNode nodei;
			final int posi;

			public Edge(KStarTreeNode nodeim1, KStarTreeNode nodei) {
				this.nodeim1 = nodeim1;
				this.posim1 = nodeim1.getAssignmentIndex();
				this.nodei = nodei;
				this.posi = nodei.getAssignmentIndex();
			}
		}

		for (KStarTreeNode topNode : collectNodesAt(root, depth - 2)) {

			// collect every edge between nodes at res i-1 and res i
			// and detatch those nodes from the tree
			List<Edge> edges = new ArrayList<>();
			for (KStarTreeNode nodeim1 : topNode.children) {
				for (KStarTreeNode nodei : nodeim1.children) {
					edges.add(new Edge(nodeim1, nodei));
				}
			}

			// remove the edges from the tree
			for (Edge edge : edges) {
				edge.nodeim1.removeFromParent();
			}

			// sanity check: the top node shouldn't have any more children left. we ate them all
			if (!topNode.children.isEmpty()) {
				throw new Error("tree levels are not fully expanded, can't sort");
			}

			// collect the nodes at res i+1 and remove them from the tree
			List<KStarTreeNode> nodesip1 = new ArrayList<>();
			for (Edge edge : edges) {
				nodesip1.addAll(edge.nodei.children);
			}
			for (KStarTreeNode nodeip1 : nodesip1) {
				nodeip1.removeFromParent();
			}

			// reverse the edge order and put the nodes back in the tree
			for (Edge edge : edges) {

				// make a new node at res i-1
				// use the assingments from the old node at res i
				// but reuse an exsting node if possible
				int posi = edge.posi;
				int rci = edge.nodei.getConfAssignments()[posi];
				String assignmenti = edge.nodei.getAssignments()[posi];
				KStarTreeNode newNodeim1 = topNode.children.stream()
					.filter(node -> node.getConfAssignments()[posi] == rci)
					.findFirst()
					.orElseGet(() -> topNode.assign(
						posi, rci, assignmenti,
						// don't add bound info yet, we'll push up from the bottom later
						BigDecimal.ZERO, BigDecimal.ZERO,
						Double.NaN, Double.NaN
					));

				// make a new node at res i
				// use the assignments from the old node at res i-1
				// use the bounds from the old node at res i
				int posim1 = edge.posim1;
				int rcim1 = edge.nodeim1.getConfAssignments()[posim1];
				String assignmentim1 = edge.nodeim1.getAssignments()[posim1];
				KStarTreeNode newNodei = newNodeim1.assign(
					posim1, rcim1, assignmentim1,
					edge.nodei.getLowerBound(), edge.nodei.getUpperBound(),
					edge.nodei.getConfLowerBound(), edge.nodei.getConfUpperBound()
				);

				// transplant the old nodes at res i+1 to be under the new node at res i
				Iterator<KStarTreeNode> iter = nodesip1.iterator();
				while (iter.hasNext()) {
					KStarTreeNode nodeip1 = iter.next();

					// does this node belong under this edge?
					int[] conf = nodeip1.getConfAssignments();
					if (conf[posim1] == rcim1 && conf[posi] == rci) {

						// yup, add it to the tree
						iter.remove();
						newNodei.addChild(nodeip1);
					}
				}
			}

			// push up bound info from nodes at res i to nodes at res i-1
			for (KStarTreeNode nodeim1 : topNode.children) {
				nodeim1.updateBoundsFromChildren(mathContext);
			}

			// sanity check: we should have transplanted all the nodes at res i+1, right?
			if (!nodesip1.isEmpty()) {
				throw new Error("tree layer swap orphaned some nodes. this is a bug.");
			}
		}
	}

	private static List<KStarTreeNode> collectNodesAt(KStarTreeNode root, int depth) {

		List<KStarTreeNode> nodes = new ArrayList<>();

		// make a recusion helper
		// (with the box hack, because javac is dumb)
		class Box<T> { public T f; }
		Box<Consumer<KStarTreeNode>> box = new Box<>();
		box.f = node -> {
			if (node.level == depth) {
				nodes.add(node);
			} else {
				for (KStarTreeNode child : node.children) {
					box.f.accept(child);
				}
			}
		};
		box.f.accept(root);

		return nodes;
	}
}
