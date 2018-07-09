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

package edu.duke.cs.osprey.externalMemory;

import java.util.*;
import java.util.function.Predicate;

import edu.duke.cs.tpie.serialization.SerializingDoublePriorityQueue;
import edu.duke.cs.tpie.serialization.SerializingFIFOQueue;

/**
 * A queue that can only be iterated once.
 */
public interface Queue<T> {
	
	void push(T val);
	T peek();
	void pop();
	long size();
	boolean isEmpty();
	
	default void pushAll(Collection<T> vals) {
		for (T val : vals) {
			push(val);
		}
	}
	
	default void pushAll(T[] vals) {
		for (T val : vals) {
			push(val);
		}
	}
	
	default T poll() {
		T val = peek();
		if (val != null) {
			pop();
		}
		return val;
	}
	
	default <Q extends Queue<T>> Q filterTo(Q other, Predicate<T> pred) {
		while (!isEmpty()) {
			T val = poll();
			if (pred.test(val)) {
				other.push(val);
			}
		}
		return other;
	}

	default Iterator<T> iterator() {
		return new Iterator<T>() {

			@Override
			public boolean hasNext() {
				return peek() != null;
			}

			@Override
			public T next() {
				if (!hasNext()) {
					throw new NoSuchElementException();
				}
				return poll();
			}
		};
	}
	
	public static interface Factory<T> {
		
		Queue<T> make();
		
		public static interface FIFO<T> {
			Queue.FIFO<T> make();
		}
	}
	
	/**
	 * A special interface for FIFO queues, since they allow
	 * more space-efficient in-place filtering.
	 */
	public static interface FIFO<T> extends Queue<T> {
		
		default void filter(Predicate<T> pred) {
			long n = size();
			for (long i=0; i<n; i++) {
				T val = poll();
				if (pred.test(val)) {
					push(val);
				}
			}
		}
	}
	
	public static class FIFOFactory<T> implements Factory.FIFO<T> {
		
		@SafeVarargs
		public static <T> Queue.FIFO<T> of(T ... vals) {
			return new Queue.FIFO<T>() {
				
				private ArrayDeque<T> q;
				{
					q = new ArrayDeque<>();
					for (T val : vals) {
						q.add(val);
					}
				}

				@Override
				public void push(T val) {
					q.add(val);
				}

				@Override
				public T peek() {
					return q.peek();
				}

				@Override
				public void pop() {
					q.pop();
				}

				@Override
				public long size() {
					return q.size();
				}

				@Override
				public boolean isEmpty() {
					return q.isEmpty();
				}
			};
		}

		public static <T> Queue.FIFO<T> of(Iterable<T> iterable, long size) {
			return of(iterable.iterator(), size);
		}

		public static <T> Queue.FIFO<T> of(Iterator<T> iter, long size) {
			return new Queue.FIFO<T>() {

				private T next = null;

				{
					// get the initial element, if any
					advance();
				}

				private void advance() {
					if (iter.hasNext()) {
						next = iter.next();
					} else {
						next = null;
					}
				}

				@Override
				public void push(T val) {
					throw new UnsupportedOperationException("this FIFO queue is read only");
				}

				@Override
				public T peek() {
					return next;
				}

				@Override
				public void pop() {
					advance();
				}

				@Override
				public long size() {
					return size;
				}

				@Override
				public boolean isEmpty() {
					return next != null;
				}
			};
		}

		@Override
		public Queue.FIFO<T> make() {
			return of();
		}
	}
	
	public static class PriorityFactory<T> implements Factory<T> {
		
		@SafeVarargs
		public static <T> Queue<T> of(Comparator<T> comparator, T ... vals) {
			return new Queue<T>() {
				
				private PriorityQueue<T> q;
				{
					q = new PriorityQueue<>(comparator);
					for (T val : vals) {
						q.add(val);
					}
				}

				@Override
				public void push(T val) {
					q.add(val);
				}

				@Override
				public T peek() {
					return q.peek();
				}

				@Override
				public void pop() {
					q.remove();
					
				}

				@Override
				public long size() {
					return q.size();
				}

				@Override
				public boolean isEmpty() {
					return q.isEmpty();
				}
			};
		}
		
		private Comparator<T> comparator;
		
		public PriorityFactory(Comparator<T> comparator) {
			this.comparator = comparator;
		}

		@Override
		public Queue<T> make() {
			return of(comparator);
		}
	}
	
	public static class ExternalFIFOFactory<T> implements Factory.FIFO<T> {
		
		@SafeVarargs
		public static <T> Queue.FIFO<T> of(SerializingFIFOQueue.Serializer<T> serializer, T ... vals) {
			return new Queue.FIFO<T>() {
				
				private SerializingFIFOQueue<T> q;
				{
					q = new SerializingFIFOQueue<>(serializer);
					for (T val : vals) {
						q.push(val);
					}
				}

				@Override
				public void push(T val) {
					q.push(val);
				}

				@Override
				public T peek() {
					if (q.empty()) {
						return null;
					}
					return q.front();
				}

				@Override
				public void pop() {
					q.pop();
				}

				@Override
				public long size() {
					return q.size();
				}

				@Override
				public boolean isEmpty() {
					return q.empty();
				}
			};
		}
		
		private SerializingFIFOQueue.Serializer<T> serializer;
		
		public ExternalFIFOFactory(SerializingFIFOQueue.Serializer<T> serializer) {
			this.serializer = serializer;
		}
		
		@Override
		public Queue.FIFO<T> make() {
			return of(serializer);
		}
	}
	
	public static class ExternalPriorityFactory<T> implements Factory<T> {
		
		@SafeVarargs
		public static <T> Queue<T> of(SerializingDoublePriorityQueue.Serializer<T> serializer, T ... vals) {
			return new Queue<T>() {
				
				private SerializingDoublePriorityQueue<T> q;
				{
					q = new SerializingDoublePriorityQueue<>(serializer);
					for (T val : vals) {
						q.push(val);
					}
				}

				@Override
				public void push(T val) {
					q.push(val);
				}

				@Override
				public T peek() {
					return q.top();
				}

				@Override
				public void pop() {
					q.pop();
				}

				@Override
				public long size() {
					return q.size();
				}

				@Override
				public boolean isEmpty() {
					return q.empty();
				}
			};
		}
		
		private SerializingDoublePriorityQueue.Serializer<T> serializer;
		
		public ExternalPriorityFactory(SerializingDoublePriorityQueue.Serializer<T> serializer) {
			this.serializer = serializer;
		}
		
		@Override
		public Queue<T> make() {
			return of(serializer);
		}
	}
}
