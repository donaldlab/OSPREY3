package edu.duke.cs.osprey.tools;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.Test;


public class TestJvmMem {

	@Test
	public void test() {

		assertThat(JvmMem.getEdenPool().name, endsWith("Eden Space"));
		assertThat(JvmMem.getSurvivorPool().name, endsWith("Survivor Space"));
		assertThat(JvmMem.getOldPool().name, endsWith("Old Gen"));
	}
}
