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

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import cern.jet.math.Constants;

public class Profiler {
	
	private Map<String,Stopwatch> stopwatches;
	private Stopwatch current;
	
	public Profiler() {
		stopwatches = new LinkedHashMap<>();
		current = null;
	}
	
	public Profiler(String name) {
		this();
		start(name);
	}
	
	public Stopwatch getOrMake(String name) {
		Stopwatch stopwatch = stopwatches.get(name);
		if (stopwatch == null) {
			stopwatch = new Stopwatch();
			stopwatches.put(name, stopwatch);
		}
		return stopwatch;
	}
	
	public void start(String name) {
		stop();
		current = getOrMake(name);
		current.start();
	}
	
	public void resume(String name) {
		stop();
		current = getOrMake(name);
		current.resume();
	}
	
	public void stop() {
		if (current != null) {
			current.stop();
			current = null;
		}
	}
	
	public String makeReport(TimeUnit timeUnit) {
		
		stop();
		
		StringBuilder buf = new StringBuilder();
		buf.append("profiling:");
		long totalNs = 0;
		for (Map.Entry<String,Stopwatch> entry : stopwatches.entrySet()) {
			buf.append(String.format("   %s: %8s", entry.getKey(), entry.getValue().getTime(timeUnit)));
			totalNs += entry.getValue().getTimeNs();
		}
		buf.append(String.format("   %s: %8s", "total", TimeFormatter.format(totalNs, timeUnit)));
		buf.append(String.format(" %10.2f ops", (float)TimeFormatter.NSpS/totalNs));
		return buf.toString();
	}
}
