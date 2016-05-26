package edu.duke.cs.osprey.tools;

import java.util.concurrent.TimeUnit;

public class Stopwatch {

	public static final long NSpUS = 1000;
	public static final long NSpMS = NSpUS * 1000;
	public static final long NSpS = NSpMS * 1000;
	public static final long NSpM = NSpS * 60;
	public static final long NSpH = NSpM * 60;
	
	private boolean isRunning;
	private long startTime;
	private long stopTime;
	
	public Stopwatch() {
		isRunning = false;
		startTime = -1;
	}
	
	public void start() {
		assert (!isRunning);
		startTime = System.nanoTime();
		isRunning = true;
	}
	
	public void stop() {
		assert (isRunning);
		stopTime = System.nanoTime();
		isRunning = false;
	}
	
	public boolean isRunning() {
		return isRunning;
	}
	
	public long getTimeNs() {
		if (isRunning) {
			return System.nanoTime() - startTime;
		} else {
			return stopTime - startTime;
		}
	}
	
	public double getTimeUs() {
		return (double)getTimeNs()/NSpUS;
	}
	
	public double getTimeMs() {
		return (double)getTimeNs()/NSpMS;
	}
	
	public double getTimeS() {
		return (double)getTimeNs()/NSpS;
	}
	
	public double getTimeM() {
		return (double)getTimeNs()/NSpM;
	}
	
	public double getTimeH() {
		return (double)getTimeNs()/NSpH;
	}
	
	public String getTime() {
		return getTime(0);
	}
	
	public String getTime(int decimals) {
		if (getTimeH() > 1) {
			return getTime(TimeUnit.HOURS, decimals);
		} else if (getTimeM() > 1) {
			return getTime(TimeUnit.MINUTES, decimals);
		} else if (getTimeS() > 1) {
			return getTime(TimeUnit.SECONDS, decimals);
		} else if (getTimeMs() > 1) {
			return getTime(TimeUnit.MILLISECONDS, decimals);
		} else if (getTimeUs() > 1) {
			return getTime(TimeUnit.MICROSECONDS, decimals);
		} else {
			return getTime(TimeUnit.NANOSECONDS, decimals);
		}
	}
	
	public String getTime(TimeUnit unit) {
		return getTime(unit, 0);
	}
	
	public String getTime(TimeUnit unit, int decimals) {
		String formatSpec = "%." + decimals + "f ";
		switch (unit) {
			case NANOSECONDS: return String.format(formatSpec + "ns", getTimeNs());
			case MICROSECONDS: return String.format(formatSpec + "us", getTimeUs());
			case MILLISECONDS: return String.format(formatSpec + "ms", getTimeMs());
			case SECONDS: return String.format(formatSpec + "s", getTimeS());
			case MINUTES: return String.format(formatSpec + "m", getTimeM());
			case HOURS: return String.format(formatSpec + "h", getTimeH());
			default:
				throw new Error("Unsupported time unit: " + unit);
		}
	}
}
