package edu.duke.cs.osprey.tools;

public class TimeTools {

	private static final long baseMs = System.currentTimeMillis();
	private static final long baseNs = System.nanoTime();

	/**
	 * Returns the number of milliseconds elapsed since midnight, January 1, 1970 UTC.
	 *
	 * Won't overflow a signed long until about the year 292473178. Hopefully by then we'll still a have a sun.
	 */
	public static long getTimestampMs() {
		return System.currentTimeMillis();
	}

	/**
	 * Returns approximately the number of microseconds elapsed since midnight, January 1, 1970 UTC.
	 *
	 * Won't overflow a signed long until about the year 294441. Hopefully by then we won't be using this code anymore.
	 */
	public static long getTimestampUs() {
		return baseMs*1000L + (System.nanoTime() - baseNs)/1000;
	}

	/**
	 * Returns approximately the number of nanoseconds elapsed since midnight, January 1, 1970 UTC.
	 *
	 * Won't overflow a signed long until about the year 2262. Hopefully by then we'll have better computers.
	 * And flying cars.
	 */
	public static long getTimestampNs() {
		return baseMs*1000000L + System.nanoTime() - baseNs;
	}
}
