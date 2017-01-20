package edu.duke.cs.osprey;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

public class ScriptPlayground {
	
	public static void main(String[] args)
	throws Exception {
		minimalGMEC();
	}
	
	private static void minimalGMEC()
	throws Exception {
		
		Molecule mol = PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb");
		Strand strand = new Strand(mol, "2", "73");
	}
}
