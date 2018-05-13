/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
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
