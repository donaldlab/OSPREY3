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

package edu.duke.cs.osprey.lute;


import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.Test;

import java.util.Arrays;


public class TestRandomizedDFSSampler {

	/*
	 * Trying to trigger this exception, which was observed in the wild:
	 *
	 * edu.duke.cs.osprey.confspace.TuplesIndex$NoSuchTupleException: no tuple found matching [3=0,2=0]
	 *   at edu.duke.cs.osprey.confspace.TuplesIndex.forEachIn(TuplesIndex.java:189)
	 *   at edu.duke.cs.osprey.lute.ConfSampler$Samples.addConf(ConfSampler.java:137)
	 *   at edu.duke.cs.osprey.lute.RandomizedDFSConfSampler.sampleConfsForTuples(RandomizedDFSConfSampler.java:57)
	 *   at edu.duke.cs.osprey.lute.LUTE.fit(LUTE.java:733)
	 *   at edu.duke.cs.osprey.lute.LUTE.sampleTuplesAndFit(LUTE.java:701)
	 *
	 * So far, no luck...
	 */
	@Test(expected = ConfSampler.NoMoreSamplesException.class) // we're supposed to run out of tuples to sample
	public void test() {

		// configure the forcefield
		ForcefieldParams ffparams = new ForcefieldParams();

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

		// make sure all strands share the same template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("G648", "G654")
			.build();
		for (String resNum : Arrays.asList("G648", "G649", "G650")) {
			protein.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		}

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		for (String resNum : Arrays.asList("A155", "A156", "A157")) {
			ligand.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		}

		// make the conf space
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrands(protein, ligand)
			.build();

		// compute the energy matrix
		EnergyMatrix emat;
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

			emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.build()
				.calcEnergyMatrix();
		}

		// compute the pruning matrix
		PruningMatrix pmat = new SimpleDEE.Runner()
			.setSinglesThreshold(100.0)
			.setPairsThreshold(100.0)
			.setGoldsteinDiffThreshold(20.0)
			.setTransitivePruning(true)
			.run(confSpace, emat);

		// gather some tuples
		LUTE lute = new LUTE(confSpace);
		lute.addTuples(lute.getUnprunedPairTuples(pmat));
		for (int i=0; i<3; i++) {
			lute.addUniqueTuples(lute.sampleTripleTuplesByStrongInteractions(emat, pmat, i));
		}

		// sample tuples until we can't anymore
		final int randomSeed = 12345;
		ConfSampler sampler = new RandomizedDFSConfSampler(confSpace, pmat, randomSeed);
		for (int i=1; i<Integer.MAX_VALUE; i++) {
			sampler.sampleConfsForTuples(lute.trainingSet, i++);
		}
	}
}
