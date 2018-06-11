package edu.duke.cs.osprey.newEwakstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.Comets;
import edu.duke.cs.osprey.lute.*;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import static edu.duke.cs.osprey.tools.Log.log;


public class EwakstarLab {

	public static void main(String[] args) {

		// configure the forcefield
		ForcefieldParams ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

		// make sure all strands share the same template library
		// especially since all the state conf spaces will add/share wild-type rotamers to/from the library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("G648", "G654")
				.build();
		protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "ILE", "VAL", "LEU", "ALA", "TYR").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A155", "A194")
				.build();
		ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the COMETS states
		Ewakstar.State P = new Ewakstar.State(
				"P",
				new SimpleConfSpace.Builder()
						.addStrand(protein)
						.build()
		);

		Ewakstar.State L = new Ewakstar.State(
				"L",
				new SimpleConfSpace.Builder()
						.addStrand(ligand)
						.build()
		);

		Ewakstar.State PL = new Ewakstar.State(
				"PL",
				new SimpleConfSpace.Builder()
						.addStrand(protein)
						.addStrand(ligand)
						.build()
		);


		Ewakstar ewakstar = new Ewakstar.Builder()
				.addState(PL)
				.setBoundEw(30)
				.setMaxSimultaneousMutations(1)
				.setLogFile(new File("ewakstar.sequences.tsv"))
				.build();

		// make the ecalc from all the conf spaces
		List<SimpleConfSpace> confSpaces = Arrays.asList(P.confSpace, L.confSpace, PL.confSpace);

		EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces, ffparams)
				.setParallelism(Parallelism.makeCpu(4))
				.build();

		// what are conformation energies?
		SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(PL.confSpace, ecalc)
				.build()
				.calcReferenceEnergies();
		PL.confEcalc = new ConfEnergyCalculator.Builder(PL.confSpace, ecalc)
				.setReferenceEnergies(eref)
				.build();

		// calc the energy matrix
		EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(PL.confEcalc)
				.setCacheFile(new File(String.format("emat.%s.dat", "PL")))
				.build()
				.calcEnergyMatrix();
		PL.fragmentEnergies = emat;

		// run DEE (super important for good LUTE fits!!)
		PruningMatrix pmat = new SimpleDEE.Runner()
				.setGoldsteinDiffThreshold(10.0)
				.setTypeDependent(true)
				.setShowProgress(true)
				.setCacheFile(new File(String.format("ewakstar.%s.pmat.dat", PL.name)))
				.setParallelism(Parallelism.makeCpu(8))
				.run(PL.confSpace, emat);


		// make the conf tree factory
		PL.confTreeFactory = (rcs) -> new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build();

		// use ConfDBs
		PL.confDBFile = new File("conf." + PL.name.toLowerCase() + ".db");


		// run COMETS
		List<Ewakstar.SequenceInfo> seqs = ewakstar.findBestSequences(6, 10, PL);

		for (Ewakstar.SequenceInfo info : seqs) {
			log("sequence:   %s", info.sequence);
		}
	}

}
