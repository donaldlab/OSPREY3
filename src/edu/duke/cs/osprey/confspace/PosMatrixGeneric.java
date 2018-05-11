package edu.duke.cs.osprey.confspace;

public class PosMatrixGeneric<T> extends AbstractPosMatrix<T> {

	private T[] vals;

	public PosMatrixGeneric(SimpleConfSpace confSpace) {
		super(confSpace);
	}

	public PosMatrixGeneric(int numPos) {
		super(numPos);
	}

	@Override
	@SuppressWarnings("unchecked")
	protected void allocate(int numPairs) {
		vals = (T[])new Object[numPairs];
	}

	@Override
	public T get(int pos1, int pos2) {
		return vals[getIndex(pos1, pos2)];
	}

	@Override
	public void set(int pos1, int pos2, T val) {
		vals[getIndex(pos1, pos2)] = val;
	}
}
