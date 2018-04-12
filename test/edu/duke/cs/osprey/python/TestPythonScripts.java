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
	
	private void runGMEC(String script) {
		run("examples/python.GMEC", script);
	}

	@Test public void analyzeConf() { runGMEC("analyzeConf.py"); }
	@Test public void findGMECAdvanced() { runGMEC("findGMEC.advanced.py"); }
	@Test public void findGMECDEE() { runGMEC("findGMEC.DEE.py"); }
	@Test public void findGMECDEEPer() { runGMEC("findGMEC.DEEPer.py"); }
	@Test public void findGMECExternalMemory() { runGMEC("findGMEC.externalMemory.py"); }
	@Test public void findGMECEnergyPartitions() { runGMEC("findGMEC.energyPartitions.py"); }
	@Test public void findGMECLUTE() { runGMEC("findGMEC.LUTE.py"); }
	@Test public void findGMECMultipleStrands() { runGMEC("findGMEC.multipleStrands.py"); }
	@Test public void findGMEC() { runGMEC("findGMEC.py"); }
	@Test public void findGMECReferenceEnergies() { runGMEC("findGMEC.referenceEnergies.py"); }
	@Test public void findGMECResEntropy() { runGMEC("findGMEC.resEntropy.py"); }
	@Test public void templateLibrary() { runGMEC("templateLibrary.py"); }
	@Test public void findGMECConfDB() { runGMEC("findGMEC.confDB.py"); }

	private void runKStar(String script) {
		run("examples/python.KStar", script);
	}

	@Test public void analyzeSequence() { runKStar("analyzeSequence.py"); }
	@Test public void bbkstar() { runKStar("bbkstar.py"); }
	@Test public void kstar() { runKStar("kstar.py"); }
	@Test public void kstarConfDB() { runKStar("kstar.confdb.py"); }
	@Test public void bbkstarConfDB() { runKStar("bbkstar.confdb.py"); }
}
