package edu.duke.cs.osprey.externalMemory;

import edu.duke.cs.tpie.TPIE;

public class ExternalMemory {
	
	private static boolean limitSet = false;
	
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
			throw new Error("internal memory limit already set");
		}
		TPIE.start(mib);
		limitSet = true;
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
	 * Set temporary directory to host external memory.
	 * 
	 * @param dir Path to directory. This directory must exist.
	 */
	public static void setTempDir(String dir) {
		TPIE.setTempDir(dir);
	}
	
	/**
	 * Set temporary directory to host external memory.
	 * 
	 * @param dir Path to directory. This directory must exist.
	 * @param subdir Name of subdirectory within dir. This directory will be created if it does not exist.
	 */
	public static void setTempDir(String dir, String subdir) {
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
}
