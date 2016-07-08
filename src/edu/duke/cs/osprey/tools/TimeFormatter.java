package edu.duke.cs.osprey.tools;

import java.util.concurrent.TimeUnit;

public class TimeFormatter {
	
	public static final long NSpUS = 1000;
	public static final long NSpMS = NSpUS * 1000;
	public static final long NSpS = NSpMS * 1000;
	public static final long NSpM = NSpS * 60;
	public static final long NSpH = NSpM * 60;

	public static double getTimeUs(long ns) {
		return (double)ns/NSpUS;
	}
	
	public static double getTimeMs(long ns) {
		return (double)ns/NSpMS;
	}
	
	public static double getTimeS(long ns) {
		return (double)ns/NSpS;
	}
	
	public static double getTimeM(long ns) {
		return (double)ns/NSpM;
	}
	
	public static double getTimeH(long ns) {
		return (double)ns/NSpH;
	}
	
	public static String format(long ns) {
		return format(ns, 0);
	}
	
	public static String format(long ns, int decimals) {
		long nsabs = Math.abs(ns);
		if (getTimeH(nsabs) > 1) {
			return format(ns, TimeUnit.HOURS, decimals);
		} else if (getTimeM(nsabs) > 1) {
			return format(ns, TimeUnit.MINUTES, decimals);
		} else if (getTimeS(nsabs) > 1) {
			return format(ns, TimeUnit.SECONDS, decimals);
		} else if (getTimeMs(nsabs) > 1) {
			return format(ns, TimeUnit.MILLISECONDS, decimals);
		} else if (getTimeUs(nsabs) > 1) {
			return format(ns, TimeUnit.MICROSECONDS, decimals);
		} else {
			return format(ns, TimeUnit.NANOSECONDS, decimals);
		}
	}
	
	public static String format(long ns, TimeUnit unit) {
		return format(ns, unit, 0);
	}
	
	public static String format(long ns, TimeUnit unit, int decimals) {
		
		String formatSpec;
		if (unit == TimeUnit.NANOSECONDS) {
			formatSpec = "%d ";
		} else {
			formatSpec = "%." + decimals + "f ";
		}
		
		switch (unit) {
			case NANOSECONDS: return String.format(formatSpec + "ns", ns);
			case MICROSECONDS: return String.format(formatSpec + "us", getTimeUs(ns));
			case MILLISECONDS: return String.format(formatSpec + "ms", getTimeMs(ns));
			case SECONDS: return String.format(formatSpec + "s", getTimeS(ns));
			case MINUTES: return String.format(formatSpec + "m", getTimeM(ns));
			case HOURS: return String.format(formatSpec + "h", getTimeH(ns));
			default:
				throw new Error("Unsupported time unit: " + unit);
		}
	}
}
