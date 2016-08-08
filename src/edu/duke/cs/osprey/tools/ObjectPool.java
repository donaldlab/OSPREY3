package edu.duke.cs.osprey.tools;

import java.util.ArrayDeque;
import java.util.Deque;

public class ObjectPool<T> {
	
	private Factory<T,Void> factory;
	private Deque<T> objects;
	
	public ObjectPool(Factory<T,Void> factory) {
		this.factory = factory;
		this.objects = new ArrayDeque<>();
	}
	
	public T checkout() {
		if (objects.isEmpty()) {
			objects.addLast(factory.make(null));
		}
		return objects.removeFirst();
	}
	
	public void release(T obj) {
		objects.addLast(obj);
	}
	
	public int size() {
		return objects.size();
	}
	
	public void clear() {
		objects.clear();
	}
}
