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

package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;

public class EnergyPartitionsPlayground {
	
	public static void main(String[] args) {
		
		// read a protein
		Strand strand = new Strand.Builder(PDBIO.readFile("examples/python.GMEC/1CC8.ss.pdb")).build();
		
		// configure flexibility
		for (int i=2; i<=6; i++) {
			strand.flexibility.get("A" + i).setLibraryRotamers("LEU", "ILE").setContinuous();
		}
		/*
		for (int i=2; i<=8; i++) {
			strand.flexibility.get("A" + i).setLibraryRotamers("GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "CYS", "ASN", "GLN", "ASP", "GLU", "PHE", "TRP", "TYR", "HIE", "HID", "LYS", "ARG", "MET").setContinuous();
		}
		*/
		
		// make the conf space
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand)
			//.setShellDistance(4)
			//.setShellDistance(9)
			.build();
		
		System.out.println("Residues: " + strand.mol.residues.size());
		System.out.println("Shell residues: " + confSpace.shellResNumbers.size());
		System.out.println("positions: " + confSpace.positions.size());
		System.out.println("RCs per pos: " + confSpace.positions.get(0).resConfs.size());
		
		// choose the default forcefield
		ForcefieldParams ffparams = new ForcefieldParams();
		
		EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.make(4, 1, 16))
			.setType(EnergyCalculator.Type.ResidueCudaCCD)
			.build();
		
		// get reference energies
		SimpleReferenceEnergies eref = new SimpleReferenceEnergies.Builder(confSpace, ecalc).build();
		
		// what's the energy of a conformation?
		ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
			//.setEnergyPartition(EnergyPartition.Traditional)
			.setEnergyPartition(EnergyPartition.AllOnPairs)
			.setReferenceEnergies(eref)
			.build();

		// compute the energy matrix
		EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
			.build()
			.calcEnergyMatrix();
		
		// config TPIE
		ExternalMemory.setInternalLimit(128);
		ExternalMemory.setTempDir(System.getProperty("user.home"), "tpie-epart");
		
		// how should confs be ordered?
		ConfSearch confSearch = new ConfAStarTree.Builder(emat, confSpace)
			//.setTraditional()
			.setMPLP(new ConfAStarTree.MPLPBuilder().setUpdater(new EdgeUpdater()).setNumIterations(5))
			.useExternalMemory()
			.setShowProgress(true)
			.build();
	
		// find the GMEC!
		System.out.println("Finding GMEC...");
		EnergiedConf gmec = new SimpleGMECFinder.Builder(confSearch, confEcalc)
			.setPrintIntermediateConfsToConsole(false)
			.useExternalMemory()
			.build()
			.find();
		
		/* for the original design
		// check the GMEC is still correct
		//assertThat(gmec.getAssignments(), is(new int[] { 6, 3, 9, 4, 1 }));
		//assertThat(gmec.getEnergy(), isAbsolutely(-20.891946, 1e-6));
		
		// GMEC score:      -21.623094 (gap: ~0.73)
		// min-bound score: -22.305976 (gap: ~8.89) !!!
		*/
	}
}
