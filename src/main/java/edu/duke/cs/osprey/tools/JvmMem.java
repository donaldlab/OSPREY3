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

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.lang.management.MemoryUsage;
import java.util.NoSuchElementException;

public class JvmMem {

	public static class MemInfo {

		public final String name;
		public final long usedBytes;
		public final long maxBytes;
		public final double usedPercent;

		public MemInfo(String name, long usedBytes, long maxBytes) {
			this.name = name;
			this.usedBytes = usedBytes;
			this.maxBytes = maxBytes;
			this.usedPercent = 100.0*usedBytes/maxBytes;
		}

		public MemInfo(MemoryPoolMXBean pool) {
			this(pool.getName(), pool.getUsage());
		}

		public MemInfo(String name, MemoryUsage usage) {
			this(name, usage.getUsed(), usage.getMax());
		}

		public MemInfo(MemInfo a, MemInfo b) {
			this(
				String.format("%s, %s", a.name, b.name),
				a.usedBytes + b.usedBytes,
				a.maxBytes + b.maxBytes
			);
		}

		@Override
		public String toString() {
			return String.format("%.1f%% of %s", usedPercent, MathTools.formatBytes(maxBytes));
		}
	}

	private static MemoryPoolMXBean getPool(String name) {
		return ManagementFactory.getMemoryPoolMXBeans().stream()
			.filter((pool) -> pool.getName().endsWith(name))
			.findFirst()
			.orElseThrow(() -> new NoSuchElementException("no memory pool named " + name));
	}

	/*
		On v8 JVMs with the default garbage collector (Parallel GC), the memory pool names are:
			Code Cache
			Metaspace
			Compressed Class Space
			PS Eden Space
			PS Survivor Space
			PS Old Gen

		On v11 JVMs with the default garbage collector (G1 GC), the memory pool names are:
			CodeHeap 'non-nmethods'
			Metaspace
			CodeHeap 'profiled nmethods'
			Compressed Class Space
			G1 Eden Space
			G1 Old Gen
			G1 Survivor Space
			CodeHeap 'non-profiled nmethods'
	*/

	public static MemInfo getEdenPool() {
		return new MemInfo(getPool("Eden Space"));
	}

	public static MemInfo getSurvivorPool() {
		return new MemInfo(getPool("Survivor Space"));
	}

	public static MemInfo getYoungPool() {
		return new MemInfo(getEdenPool(), getSurvivorPool());
	}

	public static MemInfo getOldPool() {
		return new MemInfo(getPool("Old Gen"));
	}
}
