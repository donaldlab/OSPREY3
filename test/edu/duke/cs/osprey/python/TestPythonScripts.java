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
	@Test public void findGMECLUTE() { runGMEC("findGMEC.LUTE.train.py"); runGMEC("findGMEC.LUTE.design.py"); }
	@Test public void findGMECMultipleStrands() { runGMEC("findGMEC.multipleStrands.py"); }
	@Test public void findGMEC() { runGMEC("findGMEC.py"); }
	@Test public void findGMECReferenceEnergies() { runGMEC("findGMEC.referenceEnergies.py"); }
	@Test public void findGMECResEntropy() { runGMEC("findGMEC.resEntropy.py"); }
	@Test public void templateLibrary() { runGMEC("templateLibrary.py"); }
	@Test public void findGMECConfDB() { runGMEC("findGMEC.confDB.py"); }
	@Test public void comets() { runGMEC("comets.py"); }
	@Test public void cometsBoundedMemory() { runGMEC("comets.boudedMemory.py"); }

	private void runKStar(String script) {
		run("examples/python.KStar", script);
	}

	@Test public void analyzeSequence() { runKStar("analyzeSequence.py"); }
	@Test public void bbkstar() { runKStar("bbkstar.py"); }
	@Test public void kstar() { runKStar("kstar.py"); }
	@Test public void kstarConfDB() { runKStar("kstar.confdb.py"); }
	@Test public void bbkstarConfDB() { runKStar("bbkstar.confdb.py"); }
	@Test public void LUTE() { runKStar("LUTE.train.py"); runKStar("LUTE.kstar.py"); runKStar("LUTE.bbkstar.py"); }
	@Test public void mskstar() { runKStar("mskstar.py"); }
	@Test public void mskstarBoudedMemory() { runKStar("mskstar.boundedMemory.py"); }

	private void runSofea(String script) {
		run("examples/python.SOFEA", script);
	}

	@Test public void sofea() { runSofea("sofea.py"); }
}
