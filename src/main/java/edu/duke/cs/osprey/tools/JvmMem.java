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
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.time.Duration;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;

public class JvmMem {

	private static boolean warnedAboutLimits = false;

	public static class MemInfo {

		public final String name;

		/** the limit, if any, on this memory pool */
		public final long maxBytes;

		/** number of bytes used by the JVM */
		public final long usedBytes;
		public final double usedPercent;

		/** number of bytes allocated by the JVM from the OS */
		public final long allocatedBytes;
		public final double allocatedPercent;

		public MemInfo(String name, long maxBytes, long usedBytes, long allocatedBytes) {
			this.name = name;
			this.maxBytes = maxBytes;
			this.usedBytes = usedBytes;
			this.usedPercent = 100.0*usedBytes/maxBytes;
			this.allocatedBytes = allocatedBytes;
			this.allocatedPercent = 100.0*allocatedBytes/maxBytes;
		}

		public MemInfo(MemoryPoolMXBean pool) {
			this(pool.getName(), List.of(pool));
		}

		public MemInfo(String name, List<MemoryPoolMXBean> pools) {
			this(
				name,
				sumMax(pools),
				pools.stream().mapToLong(it -> it.getUsage().getUsed()).sum(),
				pools.stream().mapToLong(it -> it.getUsage().getCommitted()).sum()
			);
		}

		private static long sumMax(List<MemoryPoolMXBean> pools) {

			if (pools.stream().anyMatch(it -> it.getUsage().getMax() < 0)) {
				return -1;
			}

			return pools.stream().mapToLong(it -> it.getUsage().getMax()).sum();
		}

		@Override
		public String toString() {
			if (maxBytes >= 0) {
				return String.format("%.1f%% of %s, %s alloc",
					usedPercent,
					MathTools.formatBytes(maxBytes),
					MathTools.formatBytes(allocatedBytes)
				);
			} else {
				return String.format("%s, %s alloc",
					MathTools.formatBytes(usedBytes),
					MathTools.formatBytes(allocatedBytes)
				);
			}
		}
	}

	private static List<MemoryPoolMXBean> getPools(List<String> name) {
		return ManagementFactory.getMemoryPoolMXBeans().stream()
			.filter(pool ->
				name.stream()
					.anyMatch(n -> pool.getName().endsWith(n))
			)
			.collect(Collectors.toList());
	}

	private static MemoryPoolMXBean getPool(String name) {
		return getPools(List.of(name)).stream()
			.findFirst()
			.orElseThrow(() -> new NoSuchElementException("no memory pool named " + name));
	}

	private static List<MemoryPoolMXBean> getOtherPools(List<String> name) {
		return ManagementFactory.getMemoryPoolMXBeans().stream()
			.filter(pool ->
				name.stream()
					.noneMatch(n -> pool.getName().endsWith(n))
			)
			.collect(Collectors.toList());
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

		v14 JVMs seem have the same default GC settings as the v11 JVMs
	*/

	private static final String EdenSuffix = "Eden Space";
	private static final String SurvivorSuffix = "Survivor Space";
	private static final String OldSuffix = "Old Gen";

	public static MemInfo getEdenPool() {
		return new MemInfo(getPool(EdenSuffix));
	}

	public static MemInfo getSurvivorPool() {
		return new MemInfo(getPool(SurvivorSuffix));
	}

	public static MemInfo getYoungPool() {
		return new MemInfo("Young Gen", getPools(List.of(
			EdenSuffix, SurvivorSuffix
		)));
	}

	public static MemInfo getOldPool() {
		return new MemInfo(getPool(OldSuffix));
	}

	public static MemInfo getOverheadPool() {
		return new MemInfo("Overhead", getOtherPools(List.of(
			EdenSuffix, SurvivorSuffix, OldSuffix
		)));
	}

	public static MemInfo getProcess() {

		var young = getYoungPool();
		var old = getOldPool();
		var overhead = getOverheadPool();

		// look up memory limits for this process, if any
		// eg, SLURM will set limits on RSS
		long max = -1;
		try {
			max = GLibC.getrlimit(GLibC.RLimitResource.RSS).max;
		} catch (Throwable t) {
			if (!warnedAboutLimits) {
				System.err.println("WARN: can't determine process memory limits");
				t.printStackTrace(System.err);
				warnedAboutLimits = true;
			}
		}

		return new MemInfo(
			"Total",
			max,
			young.usedBytes + old.usedBytes + overhead.usedBytes,
			young.allocatedBytes + old.allocatedBytes + overhead.allocatedBytes
		);
	}

	/**
	 * A one line string memory report that shows usage and limits for the JVM heap and the JVM process.
	 */
	public static String oneLineReport() {
		var old = getOldPool();
		var process = getProcess();
		return String.format("JvmMem[Heap: %s; Process: %s]", old, process);
	}

	/**
	 * Launches a thread that periodically reports memory usage statistics to stdout.
	 * The thread will run until the JVM exits.
	 */
	public static void monitorMemoryToStdout(Duration interval, boolean showHostname) {

		// try to get the hostname, if needed
		String hostname = null;
		if (showHostname) {
			try {
				hostname = InetAddress.getLocalHost().getHostName();
			} catch (UnknownHostException ex) {
				System.err.println("WARN: can't find hostname");
				ex.printStackTrace(System.err);
			}
		}
		final String fhostname = hostname;

		Runnable report = () -> {
			var buf = new StringBuilder();
			if (fhostname != null) {
				buf.append(fhostname);
				buf.append(": ");
			}
			buf.append(oneLineReport());
			buf.append('\n');
			System.out.print(buf);
		};

		// start the monitor thread
		var thread = new Thread(() -> {
			while (true) {

				report.run();

				try {
					// silly IDE, it's not a busy wait
					//noinspection BusyWait
					Thread.sleep(interval.toMillis());
				} catch (InterruptedException ex) {
					return;
				}
			}
		});
		thread.setName("Memory Monitor");
		thread.setDaemon(true);
		thread.start();
	}
}
