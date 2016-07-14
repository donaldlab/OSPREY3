package edu.duke.cs.osprey.tools;

import java.util.Iterator;

public abstract class InfiniteIterator<T> implements Iterator<T> {

	@Override
	public boolean hasNext() {
		return true;
	}
}
