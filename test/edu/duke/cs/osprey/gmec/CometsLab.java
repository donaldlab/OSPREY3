package edu.duke.cs.osprey.gmec;

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


public class CometsLab {

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

		// make the COMETS states
		Comets.State unbound = new Comets.State(
			"Unbound",
			new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build()
		);
		Comets.State bound = new Comets.State(
			"Bound",
			new SimpleConfSpace.Builder()
				.addStrand(protein)
				.addStrand(ligand)
				.build()
		);

		// configure COMETS
		Comets.LME objective = new Comets.LME.Builder()
			.addState(bound, 1.0)
			.addState(unbound, -1.0)
			.build();
		Comets.LME constraint = new Comets.LME.Builder()
			.addState(bound, 1.0)
			.addState(unbound, -1.0)
			.constrainLessThan(-20.0)
			.build();
		Comets comets = new Comets.Builder(objective)
			//.addConstraint(constraint)
			.setLogFile(new File("comets.sequences.tsv"))
			.build();

		// make the ecalc from all the conf spaces
		List<SimpleConfSpace> confSpaces = Arrays.asList(unbound.confSpace, bound.confSpace);
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			// configure the COMETS states for conf searching
			for (Comets.State state : comets.states) {

				// what are conformation energies?
				SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(state.confSpace, ecalc)
					.build()
					.calcReferenceEnergies();
				state.confEcalc = new ConfEnergyCalculator.Builder(state.confSpace, ecalc)
					.setReferenceEnergies(eref)
					.build();

				// calc the energy matrix
				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(state.confEcalc)
					//.setCacheFile(new File(String.format("emat.%s.dat", state.name)))
					.build()
					.calcEnergyMatrix();
				state.fragmentEnergies = emat;

				// run DEE (super important for good LUTE fits!!)
				PruningMatrix pmat = new SimpleDEE.Runner()
					.setGoldsteinDiffThreshold(10.0)
					.setTypeDependent(true)
					.setShowProgress(true)
					//.setCacheFile(new File(String.format("comets.%s.pmat.dat", state.name)))
					.setParallelism(Parallelism.makeCpu(8))
					.run(state.confSpace, emat);

				final File luteFile = new File(String.format("comets.LUTE.%s.dat", state.name));
				if (!luteFile.exists()) {

					final File confDBFile = new File(String.format("comets.%s.conf.db", state.name));
					try (ConfDB confdb = new ConfDB(state.confSpace, confDBFile)) {
						ConfDB.ConfTable confTable = confdb.new ConfTable("COMETS");

						log("\nLUTE for state %s:\n", state.name);

						final int randomSeed = 12345;
						final LUTE.Fitter fitter = LUTE.Fitter.OLSCG;
						final double maxOverfittingScore = 1.5;
						final double maxRMSE = 0.1;

						// compute LUTE fit
						LUTE lute = new LUTE(state.confSpace);
						ConfSampler sampler = new UniformConfSampler(state.confSpace, pmat, randomSeed);
						lute.sampleTuplesAndFit(state.confEcalc, emat, pmat, confTable, sampler, fitter, maxOverfittingScore, maxRMSE);
						lute.reportConfSpaceSize(pmat);
						lute.save(luteFile);
					}
				}
				LUTEConfEnergyCalculator luteEcalc = new LUTEConfEnergyCalculator(state.confSpace, LUTEIO.read(luteFile));
				state.fragmentEnergies = luteEcalc;
				state.confEcalc = luteEcalc;

				// make the conf tree factory
				state.confTreeFactory = (rcs) -> new ConfAStarTree.Builder(emat, rcs)
					.setLUTE(luteEcalc)
					//.setTraditional()
					.build();

				// use ConfDBs
				state.confDBFile = new File("conf." + state.name.toLowerCase() + ".db");
			}

			// run COMETS
			List<Comets.SequenceInfo> seqs = comets.findBestSequences(5);

			for (Comets.SequenceInfo info : seqs) {
				log("sequence: %8.3f   %s", info.objective, info.sequence);
			}
		}
	}
}
