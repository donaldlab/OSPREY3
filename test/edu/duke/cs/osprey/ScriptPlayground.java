package edu.duke.cs.osprey;

import java.io.File;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.Forcefield;
import edu.duke.cs.osprey.gmec.LoggingConfPrinter;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;
import edu.duke.cs.osprey.minimization.SimpleConfMinimizer;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

public class ScriptPlayground {
	
	public static void main(String[] args)
	throws Exception {
		minimalGMEC();
		//advancedGMEC();
	}
	
	private static void minimalGMEC()
	throws Exception {
		
		// read a protein
		Strand strand = Strand.builder(PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb")).build();
		
		// configure flexibility
		strand.flexibility.get(2).setLibraryRotamers("ALA", "GLY");
		strand.flexibility.get(3).setLibraryRotamers(Strand.WildType, "VAL");
		strand.flexibility.get(4).setLibraryRotamers();
		
		// make the conf space
		SimpleConfSpace confSpace = SimpleConfSpace.build(strand);
		
		// get an energy matrix
		EnergyMatrix emat = SimplerEnergyMatrixCalculator.build(confSpace).calcEnergyMatrix(new File("/tmp/emat.dat"));
		
		// find the GMEC!
		SimpleGMECFinder.builder(confSpace, emat).build().find(5);
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
		//protein.flexibility.get(4).setLibraryRotamers(Strand.WildType, "TRP", "PHE").setContinuousEllipses(10); // not implemented yet
		protein.flexibility.get(5).setLibraryRotamers().addWildTypeRotamers();
		
		/* TODO: make this work
		// make the ligand
		Strand ligand = Strand.builder(mol)
			.setResidues(45, 50)
			.setTemplateLibrary(templateLib)
			.build();
		
		// configure protein flexibility
		ligand.flexibility.get(45).setLibraryRotamers("PRO");
		ligand.flexibility.get(46).setLibraryRotamers(Strand.WildType, "TRP").setContinuous(4);
		*/
		
		// make the conf space
		SimpleConfSpace confSpace = SimpleConfSpace.builder()
			.addStrand(protein)
			// TEMP
			//.addStrand(ligand, new StrandFlex.TranslateRotate(10, 10))
			.build();
		
		// choose the forcefield params
		ForcefieldParams ffparams = new ForcefieldParams(ff);
		ffparams.solvationForcefield = null; // turn off solvation energy
		
		// get an energy matrix
		SimplerEnergyMatrixCalculator ematcalc = SimplerEnergyMatrixCalculator.builder(confSpace)
			.setForcefieldParams(ffparams)
			.build();
		EnergyMatrix emat = ematcalc.calcEnergyMatrix(new File("/tmp/emat.dat"));
		
		// TODO: optional pruning? or iMinDEE?
		
		// make the conf search
		ConfSearch search = ConfAStarTree.builder(emat, confSpace)
			.setMPLP(ConfAStarTree.MPLPBuilder().setNumIterations(4))
			.build();
		
		// find the GMEC!
		SimpleGMECFinder gmecFinder = SimpleGMECFinder.builder(confSpace, search)
			.setEnergyCalculator(SimpleConfMinimizer.builder(confSpace)
				.setForcefieldParams(ffparams)
				.setParallelism(Parallelism.makeCpu(2))
				.build())
			.setLogPrinter(new LoggingConfPrinter(new File("confs.txt")))
			.build();
		EnergiedConf gmec = gmecFinder.find();
		
		// save the PDB file
		PDBIO.writeFile(confSpace.makeMolecule(gmec.getAssignments()).mol, "gmec.pdb");
	}
}
