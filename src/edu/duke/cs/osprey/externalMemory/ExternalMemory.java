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

package edu.duke.cs.osprey.externalMemory;

import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.tpie.TPIE;

import java.io.File;

public class ExternalMemory {
	
	private static boolean limitSet = false;
	private static File tempDir = null;

	/**
	 * Set the maximum amount of internal memory (eg, RAM) to use for
	 * large data structures. External memory-aware data structures will
	 * use external memory (eg, disk, SSD, NAS) for extra storage space
	 * when internal memory limits have been reached.
	 *  
	 * @param mib maximum amount of internal memory to use, in MiB
	 */
	public static void setInternalLimit(int mib) {
		if (limitSet) {
			System.err.println("WARNING: Internal memory limit already set, ignoring additional request.");
			return;
		}
		TPIE.start(mib);
		limitSet = true;
		setDefaultTempDir();
	}
	
	/**
	 * Throw a {@link InternalMemoryLimitNotSetException} if the internal memory limit has not yet been set by {@link #setInternalLimit(int)}.
	 */
	public static void checkInternalLimitSet() {
		if (!limitSet) {
			throw new InternalMemoryLimitNotSetException();
		}
	}
	
	public static class InternalMemoryLimitNotSetException extends RuntimeException {

		private static final long serialVersionUID = 1672062453287951982L;
		
		public InternalMemoryLimitNotSetException() {
			super(
				"Before external-memory aware data structures can be used,"
				+ " the internal memory limit must be set."
				+ " Use ExternalMemory.setInternalLimit() to set the limit."
			);
		}
	}

	/**
	 * Set the temporary directory for external memory to the JVM default.
	 */
	public static void setDefaultTempDir() {
		setTempDir(System.getProperty("java.io.tmpdir"));
	}
	
	/**
	 * Set temporary directory to host external memory.
	 * 
	 * @param dir Path to directory. The directory will be created if it does not exist
	 */
	public static void setTempDir(String dir) {

		// create the directory if needed
		File dirFile = new File(dir);
		if (!dirFile.exists()) {
			dirFile.mkdirs();
		}

		tempDir = new File(dir);
		TPIE.setTempDir(dir);
	}
	
	/**
	 * Set temporary directory to host external memory.
	 * 
	 * @param dir Path to directory. This directory will be created if it does not exist.
	 * @param subdir Name of subdirectory within dir. This directory will be created if it does not exist.
	 */
	public static void setTempDir(String dir, String subdir) {

		// create the directory if needed
		File dirFile = new File(dir);
		if (!dirFile.exists()) {
			dirFile.mkdirs();
		}

		tempDir = new File(dir);
		TPIE.setTempDir(dir, subdir);
	}
	
	/**
	 * Return true if the internal memory limit has been set by a call to {@link #setInternalLimit(int)}.
	 * @return
	 */
	public static boolean isInternalLimitSet() {
		return limitSet;
	}
	
	/**
	 * Return the number of bytes currently used in external memory. (ie, written to the temporary directory)
	 */
	public static long getExternalBytes() {
		if (!limitSet) {
			return 0;
		}
		return TPIE.getExternalBytes();
	}

	public static String getUsageReport() {
		long usedBytes = getExternalBytes();
		if (tempDir != null) {
			long freeBytes = tempDir.getUsableSpace();
			return String.format("%s (%.1f%% of %s)",
				MathTools.formatBytes(usedBytes),
				100.0*usedBytes/freeBytes,
				MathTools.formatBytes(freeBytes)
			);
		} else {
			return MathTools.formatBytes(usedBytes);
		}
	}
	
	/**
	 * Clean up resources used by the external memory system.
	 * 
	 * Under normal circumstances, this will be called automatically when the JVM exits,
	 * and you won't have to call it manually.
	 */
	public static void cleanup() {
		TPIE.stop();
		limitSet = false;
		tempDir = null;
	}
	
	/**
	 * Convenience method to initialize the external memory system, run a block of code,
	 * and then make sure the external memory gets cleaned up before returning.
	 * 
	 * @param internalMiB maximum amount of internal memory to use, in MiB 
	 * @param block A block of code to run using external memory.
	 */
	public static void use(int internalMiB, TPIE.Block block) {
		limitSet = true;
		try {
			TPIE.use(internalMiB, () -> {
				setDefaultTempDir();
				block.run();
			});
		} finally {
			limitSet = false;
			tempDir = null;
		}
	}
}
