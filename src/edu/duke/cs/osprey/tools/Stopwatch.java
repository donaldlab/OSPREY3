package edu.duke.cs.osprey.tools;

import java.util.concurrent.TimeUnit;

public class Stopwatch {

	public static final long NSpUS = 1000;
	public static final long NSpMS = 1000 * 1000;
	public static final long NSpS = 1000 * 1000 * 1000;
	
	private static boolean isRunning;
	private static long startTime;
	private static long stopTime;
	
	static {
		isRunning = false;
		startTime = -1;
	}
	
	public static void start() {
		assert (!isRunning);
		startTime = System.nanoTime();
		isRunning = true;
	}
	
	public static void stop() {
		assert (isRunning);
		stopTime = System.nanoTime();
		isRunning = false;
	}
	
	public static long getTimeNs() {
		assert (!isRunning);
		return stopTime - startTime;
	}
	
	public static double getTimeUs() {
		return (double)getTimeNs()/NSpUS;
	}
	
	public static double getTimeMs() {
		return (double)getTimeNs()/NSpMS;
	}
	
	public static double getTimeS() {
		return (double)getTimeNs()/NSpS;
	}
	
	public static String getTime(TimeUnit unit) {
		return getTime(unit, 0);
	}
	
	public static String getTime(TimeUnit unit, int decimals) {
		String formatSpec = "%." + decimals + "f ";
		switch (unit) {
			case NANOSECONDS: return String.format(formatSpec + "ns", getTimeNs());
			case MICROSECONDS: return String.format(formatSpec + "us", getTimeUs());
			case MILLISECONDS: return String.format(formatSpec + "ms", getTimeMs());
			case SECONDS: return String.format(formatSpec + "s", getTimeS());
			default:
				throw new Error("Unsupported time unit: " + unit);
		}
	}
}
