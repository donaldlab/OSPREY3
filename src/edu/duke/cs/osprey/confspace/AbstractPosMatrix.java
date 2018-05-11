package edu.duke.cs.osprey.confspace;

import java.io.Serializable;
import java.util.Iterator;

public abstract class AbstractPosMatrix<T> implements PosMatrix<T>, Serializable {

	private int numPos;
	private int numPairs;

	protected AbstractPosMatrix(SimpleConfSpace confSpace) {
		this(confSpace.positions.size());
	}

	protected AbstractPosMatrix(int numPos) {
		this.numPos = numPos;
		this.numPairs = numPos*(numPos - 1)/2;
		allocate(numPairs);
	}

	protected abstract void allocate(int numPairs);

	@Override
	public int getNumPos() {
		return numPos;
	}

	@Override
	public int getNumPairs() {
		return numPairs;
	}

	private int getIndexNoCheck(int pos1, int pos2) {
		return pos1*(pos1 - 1)/2 + pos2;
	}

	protected int getIndex(int pos1, int pos2) {
		
		// res2 should be strictly less than res1
		if (pos2 > pos1) {
			int swap = pos1;
			pos1 = pos2;
			pos2 = swap;
		} else if (pos1 == pos2) {
			throw new Error("Can't pair residue " + pos1 + " with itself");
		}
		
		return getIndexNoCheck(pos1, pos2);
	}
	
	@Override
	public void fill(T val) {
		for (int pos1=0; pos1<getNumPos(); pos1++) {
			for (int pos2=0; pos2<pos1; pos2++) {
				set(pos1, pos2, val);
			}
		}
	}
	
	@Override
	public void fill(Iterator<T> val) {
		for (int pos1=0; pos1<getNumPos(); pos1++) {
			for (int pos2=0; pos2<pos1; pos2++) {
				set(pos1, pos2, val.next());
			}
		}
	}
}
