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

package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.util.Arrays;

public class BenchmarkGMECFinder {

	public static void main(String[] args) {

		Molecule mol = PDBIO.readFile("examples/python.GMEC/1CC8.ss.pdb");

		// make a conf space
		Strand strand = new Strand.Builder(mol).build();
		for (String resNum : Arrays.asList("A2", "A3", "A4", "A5", "A6", "A7", "A9", "A10", "A11", "A12")) {
			strand.flexibility.get(resNum).setLibraryRotamers("GLY", "ALA", "VAL", "LEU").addWildTypeRotamers().setContinuous();
		}
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

		Parallelism parallelism = Parallelism.make(6, 0, 0);

		// calc the emat
		EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(parallelism)
			.build();
		ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();
		EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
			.build()
			.calcEnergyMatrix();

		final double energyWindowSize = 10.0;

		// do some warmup
		new SimpleGMECFinder.Builder(
				new ConfAStarTree.Builder(emat, confSpace).build(),
				confEcalc
			)
			.build()
			.find(1.0);

		// benchmark vanilla finder
		Stopwatch vanillaSw = new Stopwatch().start();
		new SimpleGMECFinder.Builder(
				new ConfAStarTree.Builder(emat, confSpace).build(),
				confEcalc
			)
			.build()
			.find(energyWindowSize);
		System.out.println("\n\nVanilla: " + vanillaSw.stop().getTime(2));

		// benchmark confdb finder
		File confdbFile = new File("conf.db");
		confdbFile.delete();
		Stopwatch confdbSw = new Stopwatch().start();
		new SimpleGMECFinder.Builder(
				new ConfAStarTree.Builder(emat, confSpace).build(),
				confEcalc
			)
			.setConfDB(confdbFile)
			.build()
			.find(energyWindowSize);
		System.out.println("\n\nConfDB: " + confdbSw.stop().getTime(2));
		confdbFile.delete();
	}
}
