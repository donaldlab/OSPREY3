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

import java.util.concurrent.TimeUnit;

public class Stopwatch {

	private boolean isRunning;
	private long startTime;
	private long timeNs;
	
	public Stopwatch() {
		isRunning = false;
		startTime = -1;
	}
	
	public Stopwatch reset() {
		timeNs = 0;
		return this;
	}
	
	public Stopwatch start() {
		reset();
		return resume();
	}
	
	public Stopwatch resume() {
		assert (!isRunning);
		startTime = System.nanoTime();
		isRunning = true;
		return this;
	}
	
	public Stopwatch stop() {
		assert (isRunning);
		timeNs += System.nanoTime() - startTime;
		isRunning = false;
		return this;
	}
	
	public boolean isRunning() {
		return isRunning;
	}
	
	public long getTimeNs() {
		if (isRunning) {
			return timeNs + System.nanoTime() - startTime;
		} else {
			return timeNs;
		}
	}
	
	public double getTimeUs() {
		return TimeFormatter.getTimeUs(getTimeNs());
	}
	
	public double getTimeMs() {
		return TimeFormatter.getTimeMs(getTimeNs());
	}
	
	public double getTimeS() {
		return TimeFormatter.getTimeS(getTimeNs());
	}
	
	public double getTimeM() {
		return TimeFormatter.getTimeM(getTimeNs());
	}
	
	public double getTimeH() {
		return TimeFormatter.getTimeH(getTimeNs());
	}
	
	public String getTime() {
		return TimeFormatter.format(getTimeNs());
	}
	
	public String getTime(int decimals) {
		return TimeFormatter.format(getTimeNs(), decimals);
	}
	
	public String getTime(TimeUnit unit) {
		return TimeFormatter.format(getTimeNs(), unit);
	}
	
	public String getTime(TimeUnit unit, int decimals) {
		return TimeFormatter.format(getTimeNs(), unit, decimals);
	}
}
