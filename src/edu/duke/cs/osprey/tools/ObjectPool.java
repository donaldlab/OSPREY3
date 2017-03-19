package edu.duke.cs.osprey.tools;

import java.io.Closeable;
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
	
	public void allocate(int size) {
		while (this.size < size) {
			make();
		}
	}
	
	private void make() {
		objects.addLast(factory.make(null));
		size++;
	}
	
	public T checkout() {
		if (objects.isEmpty()) {
			make();
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
	
	public class Checkout implements Closeable {
		
		private final T thing;
		
		private Checkout() {
			synchronized(ObjectPool.this) {
				this.thing = checkout();
			}
		}

		@Override
		public void close() {
			synchronized(ObjectPool.this) {
				release(thing);
			}
		}
		
		public T get() {
			return thing;
		}
	}
	
	public Checkout autoCheckout() {
		return new Checkout();
	}
}
