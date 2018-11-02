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

package edu.duke.cs.osprey.paste;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.ReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.SVG;
import edu.duke.cs.osprey.tools.SVGPlot;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import static edu.duke.cs.osprey.tools.Log.log;


public class PasteRunner {

	public static void main(String[] args) {

		Molecule mol = PDBIO.readFile("./test-resources/4npd.A_DomainA_noH_trim_his_clean.min.pdb");
		ForcefieldParams ffparams = new ForcefieldParams();
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		Strand complex = new Strand.Builder(mol)
			.setResidues("A1", "A58")
			.setTemplateLibrary(templateLib)
			.build();
		complex.flexibility.get("A8").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		complex.flexibility.get("A41").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		complex.flexibility.get("A42").setLibraryRotamers(Strand.WildType, "THR", "LYS").addWildTypeRotamers().setContinuous();
		complex.flexibility.get("A45").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		complex.flexibility.get("A46").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();


		SimpleConfSpace proteinConfSpace = new SimpleConfSpace.Builder()
			.addStrand(complex)
			.build();

		Strand protein = new Strand.Builder(mol)
				.setResidues("A1", "A36")
				.setTemplateLibrary(templateLib)
				.build();
		protein.flexibility.get("A8").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();


		Strand ligand = new Strand.Builder(mol)
				.setResidues("A37", "A58")
				.setTemplateLibrary(templateLib)
				.build();
		ligand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A42").setLibraryRotamers(Strand.WildType, "THR", "LYS").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A46").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();


		SimpleConfSpace complexConfSpace = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.addStrand(ligand)
				.build();

		SimpleConfSpace pConfSpace = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build();

		SimpleConfSpace lConfSpace = new SimpleConfSpace.Builder()
				.addStrand(ligand)
				.build();

		runPaste(proteinConfSpace);
		//runKStar(pConfSpace, lConfSpace, complexConfSpace);
	}

	private static void runPaste(SimpleConfSpace proteinConfSpace) {

		Paste.Settings settings = new Paste.Settings.Builder()
			.setEpsilon(0.68)
			.setStabilityThreshold(null)
			.setNumMaxPfConfs(100000)
			.setUseWindowCriterion(false) //if not using window criterion make sure you have a tighter epsilon
			.setMaxSimultaneousMutations(3)
			.addScoreConsoleWriter()
			.addScoreFileWriter(new File("paste.txt"))
			.addMutFile(new File("mut.txt"))
			.build();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(proteinConfSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(1))
			.build()
		) {
			Paste paste = new Paste(proteinConfSpace, settings);

			SimpleReferenceEnergies eref = new SimpleReferenceEnergies.Builder(paste.protein.confSpace, ecalc).build();
			paste.protein.confEcalc = new ConfEnergyCalculator.Builder(paste.protein.confSpace, ecalc).setReferenceEnergies(eref).build();

			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(paste.protein.confEcalc)
					.setCacheFile(new File(String.format("paste.emat.%s.dat", paste.protein.id)))
					.build()
					.calcEnergyMatrix();

			paste.protein.confSearchFactory = (rcs) -> new ConfAStarTree.Builder(emat, rcs)
					.setTraditional()
					.build();

			paste.run();
		}
	}


	private static void runKStar(SimpleConfSpace targetConfSpace, SimpleConfSpace ligandConfSpace, SimpleConfSpace complexConfSpace) {

		KStar.Settings settings = new KStar.Settings.Builder()
				.setEpsilon(0.9999)
				.setStabilityThreshold(null)
				.setMaxSimultaneousMutations(3)
				.addScoreConsoleWriter()
				.addScoreFileWriter(new File("kstar.txt"))
				.build();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(complexConfSpace, new ForcefieldParams())
				.setParallelism(Parallelism.makeCpu(1))
				//.setParallelism(Parallelism.make(4, 1, 1))
				.build()
		) {
			KStar kstar = new KStar(targetConfSpace, ligandConfSpace, complexConfSpace, settings);
			for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

				info.confEcalc = new ConfEnergyCalculator.Builder(info.confSpace, ecalc).build();

				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
						.setCacheFile(new File(String.format("kstar.emat.%s.dat", info.id)))
						.build()
						.calcEnergyMatrix();

				info.confSearchFactory = (rcs) -> new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build();
			}
			kstar.run();
		}
	}
}
