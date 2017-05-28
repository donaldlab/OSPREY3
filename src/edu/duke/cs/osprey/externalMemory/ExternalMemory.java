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
	
	public static boolean isInternalLimitSet() {
		return limitSet;
	}
	
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
}
