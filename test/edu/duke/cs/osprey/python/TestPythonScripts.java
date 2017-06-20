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
	
	@Test
	public void findGMEC() {
		run("examples/1CC8.python", "findGMEC.py");
	}
	
	@Test
	public void findGMECAdvanced() {
		run("examples/1CC8.python", "findGMEC.advanced.py");
	}

	@Test
	public void findGMECReferenceEnergies() {
		run("examples/1CC8.python", "findGMEC.referenceEnergies.py");
	}
	
	@Test
	public void findGMECExternalMemory() {
		run("examples/1CC8.python", "findGMEC.externalMemory.py");
	}
}
