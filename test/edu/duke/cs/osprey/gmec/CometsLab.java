package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

import java.util.List;

import static edu.duke.cs.osprey.tools.Log.log;


public class CometsLab {

	public static void main(String[] args) {

		// configure the forcefield
		ForcefieldParams ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
			.addMoleculeForWildTypeRotamers(mol)
			.build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("G648", "G654")
			.build();
		//protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
		//protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "ASP").addWildTypeRotamers().setContinuous();
		//protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType, "SER", "ASN", "GLN").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		//ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		//ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		//ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the conf spaces
		SimpleConfSpace unbound = new SimpleConfSpace.Builder()
			.addStrand(protein)
			.build();
		SimpleConfSpace bound = new SimpleConfSpace.Builder()
			.addStrand(protein)
			.addStrand(ligand)
			.build();

		// configure the objective function
		Comets.LME objective = new Comets.LME.Builder()
			.addPositiveState("Bound", bound)
			.addNegativeState("Unbound", unbound)
			.build();

		// run COMETS
		Comets comets = new Comets.Builder(objective)
			// TODO: add constraints?
			.build();
		List<Comets.SequenceInfo> seqs = comets.findBestSequences(1);
	}
}
