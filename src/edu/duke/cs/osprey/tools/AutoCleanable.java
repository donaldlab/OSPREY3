package edu.duke.cs.osprey.tools;

import edu.duke.cs.tpie.Cleaner.Cleanable;

public interface AutoCleanable extends Cleanable, AutoCloseable {
	
	@Override
	default void close() {
		cleanWithoutCrashing();
	}
	
	default void cleanWithoutCrashing() {
		try {
			clean();
		} catch (Throwable t) {
			System.err.println("Exception during cleaning of " + getClass().getName() + ":");
			t.printStackTrace(System.err);
		}
	}
}
