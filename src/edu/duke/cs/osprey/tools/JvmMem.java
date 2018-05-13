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

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.lang.management.MemoryUsage;
import java.util.NoSuchElementException;

public class JvmMem {

	public static class MemInfo {

		public final long usedBytes;
		public final long maxBytes;
		public final double usedPercent;

		public MemInfo(long usedBytes, long maxBytes) {
			this.usedBytes = usedBytes;
			this.maxBytes = maxBytes;
			this.usedPercent = 100.0*usedBytes/maxBytes;
		}

		public MemInfo(MemoryPoolMXBean pool) {
			this(pool.getUsage());
		}

		public MemInfo(MemoryUsage usage) {
			this(usage.getUsed(), usage.getMax());
		}

		public MemInfo(MemoryPoolMXBean a, MemoryPoolMXBean b) {
			this(a.getUsage(), b.getUsage());
		}

		public MemInfo(MemoryUsage a, MemoryUsage b) {
			this(
				new MemInfo(a),
				new MemInfo(b)
			);
		}

		public MemInfo(MemInfo a, MemInfo b) {
			this(
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
			.filter((pool) -> pool.getName().equals(name))
			.findFirst()
			.orElseThrow(() -> new NoSuchElementException("no memory pool named " + name));
	}

	public static MemInfo getEdenPool() {
		return new MemInfo(getPool("PS Eden Space"));
	}

	public static MemInfo getSurvivorPool() {
		return new MemInfo(getPool("PS Eden Space"));
	}

	public static MemInfo getYoungPool() {
		return new MemInfo(getEdenPool(), getSurvivorPool());
	}

	public static MemInfo getOldPool() {
		return new MemInfo(getPool("PS Old Gen"));
	}
}
