/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

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
	
	public static class Checkout<T> implements Closeable {

		private ObjectPool<T> pool;
		private final T thing;
		
		private Checkout(ObjectPool<T> pool) {
			this.pool = pool;
			synchronized(pool) {
				this.thing = pool.checkout();
			}
		}

		@Override
		public void close() {
			synchronized(pool) {
				pool.release(thing);
			}
		}
		
		public T get() {
			return thing;
		}
	}
	
	public Checkout<T> autoCheckout() {
		return new Checkout<>(this);
	}
}
