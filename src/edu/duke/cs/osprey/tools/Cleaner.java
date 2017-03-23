package edu.duke.cs.osprey.tools;

import java.lang.ref.PhantomReference;
import java.lang.ref.Reference;
import java.lang.ref.ReferenceQueue;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.collections4.MultiValuedMap;
import org.apache.commons.collections4.multimap.ArrayListValuedHashMap;

public class Cleaner {
	
	public static interface Cleanable {
		void clean();
	}
	
	/*
	 * This interface just helps enforce the difference between objects that
	 * can be cleaned, and objects that allow detection of garbage collection.
	 * 
	 * Cleanable and GarbageDetectable can't be implemented by the same object!
	 * And in general, Cleanable instances should not maintain references to their enclosing
	 * GarbageDetectable.
	 * 
	 * And further than that, for objects to be GarbageDetectable, they should not
	 * leak references to other strongly-referenced objects (eg: enclosing scopes, other threads, globals).
	 * Instances of static nested classes can be leaked because the compiler doesn't create a
	 * 'this' reference to the enclosing object. Apparently lambdas that don't reference the enclosing
	 * object are ok to leak too.
	 * 
	 * Aaaaand just for more fun, the object must not have a finalizer. The finalization system
	 * apparently holds a reference to the object too.
	 * 
	 * These restrictions make the Cleaner system rather difficult to use for non-trivial objects
	 * (like thread pools). If your objects aren't getting garbage detected, chances are you've
	 * leaked a reference somewhere. You can use the heap dump option on your favorite profiler
	 * (eg: VisualVM) to find the leaked references.
	 */
	public static interface GarbageDetectable {
		// nothing to do
	}
	
	private static ReferenceQueue<GarbageDetectable> refQueue = new ReferenceQueue<>();
	private static MultiValuedMap<PhantomReference<GarbageDetectable>,Cleanable> cleanables = new ArrayListValuedHashMap<>();
	
	public static <T extends Cleanable> T addCleaner(GarbageDetectable thing, T cleanable) {
		
		// don't let 'thing' and 'cleanable' be the same
		// the phantom references won't work correctly if they are!
		// (because we keep a strong reference to cleanable)
		if (thing == cleanable) {
			throw new IllegalArgumentException("GarbageDetectable instance and Cleanable instance must be different objects");
		}
		
		// don't let the 'thing' have a finalizer
		try {
			thing.getClass().getDeclaredMethod("finalize");
			throw new IllegalArgumentException("GarbageDetectable instance cannot have a finalizer method");
		} catch (NoSuchMethodException ex) {
			// this is ok
		}
		
		synchronized (cleanables) {
			cleanables.put(new PhantomReference<>(thing, refQueue), cleanable);
		}
		
		return cleanable;
	}
	
	static {
		// start the cleaner on a background thread
		Thread thread = new Thread(() -> {
			try {
				
				List<Cleanable> refCleanables = new ArrayList<>();
				
				while (true) {
					
					// wait for a signal that an object we care about has been garbage collected
					Reference<? extends GarbageDetectable> ref = refQueue.remove();
					
					// to continue with garbage collection, we have to manually clear the reference
					ref.clear();
					
					// run any cleaners for this ref
					refCleanables.clear();
					synchronized (cleanables) {
						refCleanables.addAll(cleanables.remove(ref));
					}
					for (Cleanable refCleanable : refCleanables) {
						refCleanable.clean();
					}
					refCleanables.clear();
				}
				
			} catch (InterruptedException ex) {
				throw new Error(ex);
			}
			
		}, "Cleaner");
		thread.setDaemon(true);
		thread.start();
	}
	
	//-------------------------------------------------------
	// here's an example of how to use the cleanable system
	
	public static class ExampleResource implements Cleanable {
		
		public boolean isUsing;
		
		public ExampleResource() {
			System.out.println("Resource()   " + System.identityHashCode(this));
			isUsing = true;
		}
		
		@Override
		public void clean() {
			if (isUsing) {
				System.out.println("Resource.cleanup()   " + System.identityHashCode(this));
				isUsing = false;
			}
		}
	}
	
	public static class ExampleThing implements GarbageDetectable {
		
		private ExampleResource resource;
		
		public ExampleThing() {
			System.out.println("Thing()");
			resource = Cleaner.addCleaner(this, new ExampleResource());
		}
		
		public void cleanup() {
			System.out.println("Thing.cleanup()");
			resource.clean();
		}
	}
	
	public static void main(String[] args)
	throws Exception {
		
		// we have a thing we want to make sure gets cleaned up
		ExampleThing thing;
		
		// happy path
		thing = new ExampleThing();
		thing.cleanup();
		thing = null;
		
		// sad path
		thing = new ExampleThing();
		thing = null;
		
		// give the cleaner a chance to cleanup
		System.out.println("waiting for cleanup...");
		System.gc();
		Thread.sleep(400);
		System.out.println("done");
	}
}
