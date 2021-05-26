package edu.duke.cs.osprey.tools;

import java.util.concurrent.TimeUnit;


public class Duration {

	public final long length;
	public final TimeUnit unit;

	public Duration(long length, TimeUnit unit) {
		this.length = length;
		this.unit = unit;
	}

	public long toNanos() {
		return unit.toNanos(length);
	}
}
