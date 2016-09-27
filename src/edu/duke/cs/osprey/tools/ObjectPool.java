package edu.duke.cs.osprey.tools;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Iterator;

public class ObjectPool<T> implements Iterable<T> {
	
	private Factory<T,Void> factory;
	private Deque<T> objects;
	private int size;
	
	public ObjectPool(Factory<T,Void> factory) {
		this.factory = factory;
		this.objects = new ArrayDeque<>();
		this.size = 0;
	}
	
	public T checkout() {
		if (objects.isEmpty()) {
			objects.addLast(factory.make(null));
			size++;
		}
		return objects.removeFirst();
	}
	
	public void release(T obj) {
		objects.addLast(obj);
	}
	
	public int size() {
		return size;
	}
	
	public int available() {
		return objects.size();
	}
	
	public void clear() {
		objects.clear();
	}

	@Override
	public Iterator<T> iterator() {
		return objects.iterator();
	}
}
