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

package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunctionFactory;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

public class MSKStarLab {

	public static void main(String[] args) {

		// configure the forcefield
		ForcefieldParams ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		// define the protein strand
		Strand proteinStrand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("G648", "G654")
			.build();
		proteinStrand.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		proteinStrand.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
		proteinStrand.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "ASP").addWildTypeRotamers().setContinuous();
		proteinStrand.flexibility.get("G654").setLibraryRotamers(Strand.WildType, "SER", "ASN", "GLN").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligandStrand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		ligandStrand.flexibility.get("A156").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		ligandStrand.flexibility.get("A172").setLibraryRotamers(Strand.WildType, "ASP", "GLU").addWildTypeRotamers().setContinuous();
		ligandStrand.flexibility.get("A192").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "PHE", "TYR").addWildTypeRotamers().setContinuous();
		ligandStrand.flexibility.get("A193").setLibraryRotamers(Strand.WildType, "SER", "ASN").addWildTypeRotamers().setContinuous();

		MSKStar.State protein = new MSKStar.State(
			"Protein",
			new SimpleConfSpace.Builder()
				.addStrand(proteinStrand)
				.build()
		);
		MSKStar.State ligand = new MSKStar.State(
			"Ligand",
			new SimpleConfSpace.Builder()
				.addStrand(ligandStrand)
				.build()
		);
		MSKStar.State complex = new MSKStar.State(
			"Complex",
			new SimpleConfSpace.Builder()
				.addStrands(proteinStrand, ligandStrand)
				.build()
		);

		MSKStar.LMFE objective = new MSKStar.LMFE.Builder()
			.addState(complex, 1.0)
			.addState(protein, -1.0)
			.addState(ligand, -1.0)
			.build();
		MSKStar mskstar = new MSKStar.Builder(objective)
			.setEpsilon(0.9999)
			.setMaxSimultaneousMutations(3)
			.setLogFile(new File("mskstar.tsv"))
			.build();

		List<SimpleConfSpace> confSpaceList = mskstar.states.stream()
			.map(state -> state.confSpace)
			.collect(Collectors.toList());
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaceList, ffparams)
			.setParallelism(Parallelism.makeCpu(8))
			.build()
		) {
			for (MSKStar.State state : mskstar.states) {

				// how should we define energies of conformations?
				state.confEcalc = new ConfEnergyCalculator.Builder(state.confSpace, ecalc)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(state.confSpace, ecalc)
						.build()
						.calcReferenceEnergies()
					)
					.build();

				// calc energy matrix
				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(state.confEcalc)
					.build()
					.calcEnergyMatrix();
				state.fragmentEnergies = emat;

				// how should confs be ordered and searched?
				state.confTreeFactory = (rcs) -> new ConfAStarTree.Builder(emat, rcs)
					.setTraditional()
					.build();

				// use ConfDB
				state.confDBFile = null;// new File(String.format("mskstar.%s.conf.db", state.name.toLowerCase()));

				state.pfuncFactory = new PartitionFunctionFactory(state.confSpace, state.confEcalc, state.name);
				//state.pfuncFactory.setUseGradientDescent(state.confEcalc);
				EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpaceList, ffparams)
						.setParallelism(Parallelism.makeCpu(8))
						.setIsMinimizing(false)
						.build();
				ConfEnergyCalculator rigidConfEcalc = new ConfEnergyCalculator(state.confEcalc, rigidEcalc);
				state.pfuncFactory.setUseMARKStar(rigidConfEcalc);
			}

			mskstar.findBestSequences(5);
		}
	}
}
