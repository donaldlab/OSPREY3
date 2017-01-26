package edu.duke.cs.osprey;

import java.io.File;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.confspace.StrandFlex;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.Forcefield;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
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
		
		// read a protein
		Molecule mol = PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb");
		Strand strand = Strand.builder(mol).build();
		
		// configure flexibility
		strand.flexibility.get(2).setLibraryRotamers("ALA", "GLY");
		strand.flexibility.get(3).setLibraryRotamers();
		strand.flexibility.get(4).setLibraryRotamers();
		
		// make the conf space
		SimpleConfSpace confSpace = SimpleConfSpace.build(strand);
		
		// get an energy matrix
		SimpleEnergyMatrixCalculator ematcalc = SimpleEnergyMatrixCalculator.build(confSpace);
		EnergyMatrix emat = ematcalc.calcEnergyMatrix(new File("/tmp/emat.dat"));
	}
	
	private static void advancedGMEC()
	throws Exception {
		
		// read a PDB
		Molecule mol = PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb");
		
		// choose a forcefield
		Forcefield ff = Forcefield.AMBER;
		
		// choose a template library
		GenericResidueTemplateLibrary templateLib = GenericResidueTemplateLibrary.builder()
			.setLovellRotamers()
			.setForcefield(ff)
			.build();
		
		// make the protein
		Strand protein = Strand.builder(mol)
			.setResidues(2, 30)
			.setTemplateLibrary(templateLib)
			.build();
		
		// configure protein flexibility
		protein.flexibility.get(2).setLibraryRotamers("ALA", "GLY");
		protein.flexibility.get(3).setLibraryRotamers(Strand.WildType, "VAL", "ARG").setContinuous(10);
		protein.flexibility.get(4).setLibraryRotamers(Strand.WildType, "TRP", "PHE").setContinuousEllipses(10);
		protein.flexibility.get(5).addWildTypeRotamers();
		
		// make the ligand
		Strand ligand = Strand.builder(mol)
			.setResidues(45, 50)
			.setTemplateLibrary(templateLib)
			.build();
		
		// configure protein flexibility
		ligand.flexibility.get(45).setLibraryRotamers("PRO");
		ligand.flexibility.get(46).setLibraryRotamers(Strand.WildType, "TRP").setContinuous(4);
		
		// make the conf space
		SimpleConfSpace confSpace = SimpleConfSpace.builder()
			.addStrand(protein)
			.addStrand(ligand, new StrandFlex.TranslateRotate(10, 10))
			.build();
		
		// choose the forcefield params
		ForcefieldParams ffparams = new ForcefieldParams(ff);
		ffparams.solvationForcefield = null; // turn off solvation energy
		
		// get an energy matrix
		SimpleEnergyMatrixCalculator ematcalc = SimpleEnergyMatrixCalculator.builder(confSpace)
			.setForcefieldParams(ffparams)
			.build();
		EnergyMatrix emat = ematcalc.calcEnergyMatrix(new File("/tmp/emat.dat"));
	}
}
