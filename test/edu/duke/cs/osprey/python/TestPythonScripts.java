package edu.duke.cs.osprey.python;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.io.File;

import org.junit.Test;

public class TestPythonScripts {
	
	private void run(String dir, String script) {
		
		ProcessBuilder pb = new ProcessBuilder();
		pb.directory(new File(dir));
		pb.command("python", script);
		pb.inheritIO();
		try {
			Process p = pb.start();
			int returnCode = p.waitFor();
			assertThat("python script had an error", returnCode, is(0));
		} catch (Exception t) {
			throw new Error(t);
		}
	}
	
	private void run1CC8(String script) {
		run("examples/1CC8.python", script);
	}
	
	@Test
	public void findGMEC() {
		run1CC8("findGMEC.py");
	}
	
	@Test
	public void findGMECAdvanced() {
		run1CC8("findGMEC.advanced.py");
	}

	@Test
	public void findGMECReferenceEnergies() {
		run1CC8("findGMEC.referenceEnergies.py");
	}
	
	@Test
	public void findGMECExternalMemory() {
		run1CC8("findGMEC.externalMemory.py");
	}
	
	@Test
	public void findGMECMultipleStrands() {
		run1CC8("findGMEC.multipleStrands.py");
	}

	@Test
	public void findGMECDEEPer() {
		run1CC8("findGMEC.DEEPer.py");
	}

	@Test
	public void findGMECLUTE() {
		run1CC8("findGMEC.LUTE.py");
	}

	@Test
	public void templateLibrary() {
		run1CC8("templateLibrary.py");
	}
}
