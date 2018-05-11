package edu.duke.cs.osprey.confspace;

import java.util.Iterator;


public interface PosMatrix<T> {
	
	void fill(T val);
	void fill(Iterator<T> val);

	int getNumPos();
	int getNumPairs();

	T get(int pos1, int pos2);
	void set(int pos1, int pos2, T val);
}
