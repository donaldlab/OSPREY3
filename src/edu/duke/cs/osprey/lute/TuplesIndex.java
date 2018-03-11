package edu.duke.cs.osprey.lute;

import com.jogamp.common.util.IntObjectHashMap;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.util.*;


public class TuplesIndex {

	public static class NoSuchTupleException extends RuntimeException {

		public NoSuchTupleException(RCTuple tuple) {
			super("no tuple found matching " + tuple);
		}
	}

	private static final int TupleNotPresent = -1;

	public final SimpleConfSpace confSpace;

	private final List<RCTuple> tuples;

	private final int numPos;
	private final int[] numRCsAtPos;

	// lookup structure for dense pairs (like an energy matrix, very CPU cache friendly)
	private final int[] pairOffsets;
	private final int[] pairTupleIndices;

	// lookup structure for sparse triples (essentially a prefix tree aka "trie", not CPU cache friendly)
	private class SparseTuples {

		private static final int NotATuple = -1;

		private class Node {

			final int id;

			IntObjectHashMap children = null;
			int tupleIndex = TupleNotPresent;

			Node(int id) {
				this.id = id;
			}

			public Node makeChild(int pos, int rc) {

				if (children == null) {
					children = new IntObjectHashMap();
				}

				Node child = new Node(getId(pos, rc));
				children.put(child.id, child);
				return child;
			}

			public Node get(int pos, int rc) {
				if (children == null) {
					return null;
				}
				return (Node)children.get(getId(pos, rc));
			}
		}

		final int[] singleOffsets;
		final Node root = new Node(-1);

		public SparseTuples() {

			int numPos = confSpace.positions.size();
			singleOffsets = new int[numPos];
			int offset = 0;
			for (int pos=0; pos<numPos; pos++) {
				singleOffsets[pos] = offset;
				offset += numRCsAtPos[pos];
			}
		}

		int getId(int pos, int rc) {
			// combine the pos,rc pair into a single number, so hash table lookups are more efficient
			return singleOffsets[pos] + rc;
		}

		void add(RCTuple tuple, int index) {

			// just in case...
			tuple.checkSortedPositions();

			// hang tuples by their toes!
			// i.e., by greatest pos first
			// conf pos iteration usually follows pos3 < pos2 < pos1
			// pos1, the largest value, is always in the outer loop and should correspond to the top node
			Node parent = root;
			for (int i=tuple.size()-1; i>=0; i--) {

				int pos = tuple.pos.get(i);
				int rc = tuple.RCs.get(i);
				Node node = parent.get(pos, rc);
				if (node == null) {
					node = parent.makeChild(pos, rc);
				}

				if (i == 0) {
					node.tupleIndex = index;
				} else {
					parent = node;
				}
			}
		}
	}
	private SparseTuples sparseTuples = null;

	public TuplesIndex(SimpleConfSpace confSpace, Collection<RCTuple> tuplesCollection) {

		this.confSpace = confSpace;
		this.tuples = new ArrayList<>(tuplesCollection);

		// index the positions and RCs for fast lookups
		numPos = confSpace.positions.size();
		numRCsAtPos = confSpace.getNumResConfsByPos();

		// allocate space for the dense pair tuples index
		pairOffsets = new int[numPos*(numPos - 1)/2];
		int offset = 0;
		int index = 0;
		for (int pos1=0; pos1<numPos; pos1++) {
			for (int pos2=0; pos2<pos1; pos2++) {
				pairOffsets[index++] = offset;
				offset += numRCsAtPos[pos1]*numRCsAtPos[pos2];
			}
		}
		pairTupleIndices = new int[offset];
		Arrays.fill(pairTupleIndices, TupleNotPresent);

		// index the tuples
		for (int i=0; i<tuples.size(); i++) {
			RCTuple tuple = tuples.get(i);
			switch (tuple.size()) {

				case 2:
					setPairTupleIndex(
						tuple.pos.get(0), tuple.RCs.get(0),
						tuple.pos.get(1), tuple.RCs.get(1),
						i
					);
				break;

				case 3:
					if (sparseTuples == null) {
						sparseTuples = new SparseTuples();
					}
					sparseTuples.add(tuple, i);
				break;

				default:
					throw new UnsupportedOperationException("tuple index only supports pairs and triples for now");
			}
		}
	}

	public RCTuple get(int index) {
		return tuples.get(index);
	}

	public int size() {
		return tuples.size();
	}

	private int getPairIndex(int pos1, int rc1, int pos2, int rc2) {

		// pos2 should be strictly less than pos1
		if (pos2 > pos1) {
			int swap = pos1;
			pos1 = pos2;
			pos2 = swap;
			swap = rc1;
			rc1 = rc2;
			rc2 = swap;
		} else if (pos1 == pos2) {
			throw new Error("Can't pair residue " + pos1 + " with itself");
		}

		int pairIndex = pos1*(pos1 - 1)/2 + pos2;
		return pairOffsets[pairIndex] + numRCsAtPos[pos2]*rc1 + rc2;
	}

	private int getPairTupleIndex(int pos1, int rc1, int pos2, int rc2) {
		return pairTupleIndices[getPairIndex(pos1, rc1, pos2, rc2)];
	}

	private void setPairTupleIndex(int pos1, int rc1, int pos2, int rc2, int val) {
		pairTupleIndices[getPairIndex(pos1, rc1, pos2, rc2)] = val;
	}

	public static interface OnTuple {
		void onTuple(int index);
	}

	public void forEach(int[] conf, OnTuple callback) {

		// look for pairs
		for (int pos1=0; pos1<numPos; pos1++) {

			int rc1 = conf[pos1];
			if (rc1 == Conf.Unassigned) {
				continue;
			}

			for (int pos2=0; pos2<pos1; pos2++) {

				int rc2 = conf[pos2];
				if (rc2 == Conf.Unassigned) {
					continue;
				}

				int pairTupleIndex = getPairTupleIndex(pos1, rc1, pos2, rc2);
				if (pairTupleIndex == TupleNotPresent) {
					throw new NoSuchTupleException(new RCTuple(pos1, rc1, pos2, rc2));
				}

				callback.onTuple(pairTupleIndex);
			}
		}

		// look for triples
		if (sparseTuples != null) {

			for (int pos1=2; pos1<numPos; pos1++) {

				int rc1 = conf[pos1];
				if (rc1 == Conf.Unassigned) {
					continue;
				}

				// does any triple start with pos1,rc1?
				SparseTuples.Node node1 = sparseTuples.root.get(pos1, rc1);
				if (node1 == null) {
					continue;
				}

				for (int pos2=1; pos2<pos1; pos2++) {

					int rc2 = conf[pos2];
					if (rc2 == Conf.Unassigned) {
						continue;
					}

					// does any triple start with pos1,rc1;pos2,rc2?
					SparseTuples.Node node2 = node1.get(pos2, rc2);
					if (node2 == null) {
						continue;
					}

					for (int pos3=0; pos3<pos2; pos3++) {

						int rc3 = conf[pos3];
						if (rc3 == Conf.Unassigned) {
							continue;
						}

						// is this our triple?
						SparseTuples.Node node3 = node2.get(pos3, rc3);
						if (node3 == null) {
							continue;
						}

						// yup!
						callback.onTuple(node3.tupleIndex);
					}
				}
			}
		}
	}
}
