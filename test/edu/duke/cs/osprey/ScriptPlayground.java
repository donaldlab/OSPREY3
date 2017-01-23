package edu.duke.cs.osprey;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.FORCEFIELD;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

public class ScriptPlayground {
	
	public static void main(String[] args)
	throws Exception {
		minimalGMEC();
		advancedGMEC();
	}
	
	private static void minimalGMEC()
	throws Exception {
		
		Molecule mol = PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb");
		Strand strand = Strand.builder(mol).build();
	}
	
	private static void advancedGMEC()
	throws Exception {
		
		// read a PDB
		Molecule mol = PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb");
		
		// make the forcefield
		ForcefieldParams ffparams = new ForcefieldParams(FORCEFIELD.AMBER);
		ffparams.solvationForcefield = null;
		
		// make the strand
		Strand strand = Strand.builder(mol)
			.setResidues(2, 73)
			.setDefaultTemplateLibrary(ffparams)
			.build();
	}
}
