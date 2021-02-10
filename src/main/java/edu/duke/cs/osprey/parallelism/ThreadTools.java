package edu.duke.cs.osprey.parallelism;

import java.util.concurrent.TimeUnit;


public class ThreadTools {

	public static void sleep(int ms) {
		try {
			Thread.sleep(ms);
		} catch (InterruptedException ex) {
			throw new RuntimeException(ex);
		}
	}

	public static void sleep(long duration, TimeUnit timeUnit) {
		sleep((int)timeUnit.toMillis(duration));
	}
}
