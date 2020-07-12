package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.List;


/**
 * Uses simple flat array storing triples.
 * Probably faster than the TupleTrees that the TupleMatrices usually use.
 *
 * Hopefully we don't have more triples than can be counted in an int type.
 */
public class TripleMatrix<T> {
	
	public final int numPos;
	public final int[] numConfByPos;

	private final List<T> triples; // indices: pos i1,i2,i3 conf i1,i2,i3 where i1 > i2 > i3
	private final int[] offsets;

	public TripleMatrix(int numPos, int[] numConfByPos) {
		
		this.numPos = numPos;
		this.numConfByPos = numConfByPos;

		offsets = new int[numPos*(numPos - 1)*(numPos - 2)/6];
		int tripleoi = 0;
		int triplei = 0;
		for (int pos1=2; pos1<numPos; pos1++) {
			for (int pos2=1; pos2<pos1; pos2++) {
				for (int pos3=0; pos3<pos2; pos3++) {
					offsets[tripleoi++] = triplei;
					triplei += numConfByPos[pos1]*numConfByPos[pos2]*numConfByPos[pos3];
				}
			}
		}

		triples = new ArrayList<>(triplei);
		for (int i=0; i<triplei; i++) {
			triples.add(null);
		}
	}
	
	public TripleMatrix(ConfSpaceIteration confSpace) {
		this(confSpace.numPos(), confSpace.numConfsByPos());
	}

	public int size() {
		return triples.size();
	}

	public void fill(T val) {
		for (int i=0; i<triples.size(); i++) {
			triples.set(i, val);
		}
	}

	public int index(int posi1, int posi2, int posi3) {

		// need: posi3 < posi2 < posi1
		// sort them using a fixed swap chain
		if (posi2 > posi1) {
			int swap = posi1;
			posi1 = posi2;
			posi2 = swap;
		} else if (posi1 == posi2) {
			throw new Error("Can't pair design position " + posi1 + " with itself");
		}
		if (posi3 > posi2) {
			int swap = posi2;
			posi2 = posi3;
			posi3 = swap;
		} else if (posi2 == posi3) {
			throw new Error("Can't pair design position " + posi2 + " with itself");
		}
		if (posi2 > posi1) {
			int swap = posi1;
			posi1 = posi2;
			posi2 = swap;
		} else if (posi1 == posi2) {
			throw new Error("Can't pair design position " + posi1 + " with itself");
		}

		return posi1*(posi1 - 1)*(posi1 - 2)/6 + posi2*(posi2 - 1)/2 + posi3;
	}

	public int index(int posi1, int confi1, int posi2, int confi2, int posi3, int confi3) {

		// need: posi3 < posi2 < posi1
		// sort them using a fixed swap chain
		if (posi2 > posi1) {
			int swap = posi1;
			posi1 = posi2;
			posi2 = swap;
			swap = confi1;
			confi1 = confi2;
			confi2 = swap;
		} else if (posi1 == posi2) {
			throw new Error("Can't pair design position " + posi1 + " with itself");
		}
		if (posi3 > posi2) {
			int swap = posi2;
			posi2 = posi3;
			posi3 = swap;
			swap = confi2;
			confi2 = confi3;
			confi3 = swap;
		} else if (posi2 == posi3) {
			throw new Error("Can't pair design position " + posi2 + " with itself");
		}
		if (posi2 > posi1) {
			int swap = posi1;
			posi1 = posi2;
			posi2 = swap;
			swap = confi1;
			confi1 = confi2;
			confi2 = swap;
		} else if (posi1 == posi2) {
			throw new Error("Can't pair design position " + posi1 + " with itself");
		}

		return offsets[posi1*(posi1 - 1)*(posi1 - 2)/6 + posi2*(posi2 - 1)/2 + posi3]
			+ numConfByPos[posi3]*numConfByPos[posi2]*confi1
			+ numConfByPos[posi3]*confi2
			+ confi3;
	}

	public void set(int posi1, int confi1, int posi2, int confi2, int posi3, int confi3, T val) {
		triples.set(index(posi1, confi1, posi2, confi2, posi3, confi3), val);
	}

	public T get(int posi1, int confi1, int posi2, int confi2, int posi3, int confi3) {
		return triples.get(index(posi1, confi1, posi2, confi2, posi3, confi3));
	}
}
