package edu.duke.cs.osprey;

import java.io.File;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.Forcefield;
import edu.duke.cs.osprey.gmec.ConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.LoggingConfPrinter;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;
import edu.duke.cs.osprey.minimization.SimpleConfMinimizer;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

@SuppressWarnings("unused")
public class ScriptPlayground {
	
	public static void main(String[] args)
	throws Exception {
		minimalGMEC();
		//advancedGMEC();
	}
	
	private static void minimalGMEC()
	throws Exception {
		
		// read a protein
		Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb")).build();
		
		// configure flexibility
		strand.flexibility.get(2).setLibraryRotamers("ALA", "GLY");
		strand.flexibility.get(3).setLibraryRotamers(Strand.WildType, "VAL");
		strand.flexibility.get(4).setLibraryRotamers();
		
		// make the conf space
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
		// choose the default forcefield
		ForcefieldParams ffparams = new ForcefieldParams();
		
		// get an energy matrix
		EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ffparams)
			.setCacheFile(new File("/tmp/emat.dat"))
			.build()
			.calcEnergyMatrix();
		
		// how should confs be ordered?
		ConfSearch confSearch = new ConfAStarTree.Builder(emat, confSpace).build();
		
		// what's the energy of a conformation?
		ConfEnergyCalculator ecalc = new SimpleConfMinimizer.Builder(confSpace, ffparams).build();
		
		// find the GMEC!
		new SimpleGMECFinder.Builder(confSpace, confSearch, ecalc)
			.build()
			.find(5);
	}
	
	private static void advancedGMEC()
	throws Exception {
		
		// read a PDB
		Molecule mol = PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb");
		
		// choose a forcefield
		Forcefield ff = Forcefield.AMBER;
		
		// choose a template library
		GenericResidueTemplateLibrary templateLib = new GenericResidueTemplateLibrary.Builder()
			.setLovellRotamers()
			.setForcefield(ff)
			.build();
		
		// make the protein
		Strand protein = new Strand.Builder(mol)
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
		Strand ligand = new Strand.Builder(mol)
			.setResidues(45, 50)
			.setTemplateLibrary(templateLib)
			.build();
		
		// configure protein flexibility
		ligand.flexibility.get(45).setLibraryRotamers("PRO");
		ligand.flexibility.get(46).setLibraryRotamers(Strand.WildType, "TRP").setContinuous(4);
		*/
		
		// make the conf space
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(protein)
			// TEMP
			//.addStrand(ligand, new StrandFlex.TranslateRotate(10, 10))
			.build();
		
		// choose the forcefield params
		ForcefieldParams ffparams = new ForcefieldParams(ff);
		ffparams.solvationForcefield = null; // turn off solvation energy
		
		// get an energy matrix
		EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ffparams)
			.setCacheFile(new File("/tmp/emat.dat"))
			.build()
			.calcEnergyMatrix();
		
		// TODO: optional pruning? or iMinDEE?
		
		// how should confs be ordered?
		ConfSearch search = new ConfAStarTree.Builder(emat, confSpace)
			.setMPLP(ConfAStarTree.MPLPBuilder().setNumIterations(4))
			.build();
		
		// what's the energy of a conformation?
		ConfEnergyCalculator ecalc = new SimpleConfMinimizer.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.makeCpu(2))
			.build();
		
		// find the GMEC!
		SimpleGMECFinder gmecFinder = new SimpleGMECFinder.Builder(confSpace, search, ecalc)
			.setLogPrinter(new LoggingConfPrinter(new File("confs.txt")))
			.build();
		EnergiedConf gmec = gmecFinder.find();
		
		// save the PDB file
		PDBIO.writeFile(confSpace.makeMolecule(gmec.getAssignments()).mol, "gmec.pdb");
	}
}
