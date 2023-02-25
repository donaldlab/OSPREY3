package edu.duke.cs.osprey.tools;

import org.junit.jupiter.api.Test;


/**
 * Tests that some memory reporting fuctions run without crashing.
 */
public class TestJvmMem {

	@Test
	public void pools() {
		JvmMem.getEdenPool();
		JvmMem.getSurvivorPool();
		JvmMem.getYoungPool();
		JvmMem.getOldPool();
		JvmMem.getOverheadPool();
		JvmMem.getProcess();
	}

	@Test
	public void oneLineReport() {
		JvmMem.oneLineReport();
	}
}
