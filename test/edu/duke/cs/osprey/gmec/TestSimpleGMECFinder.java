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

import static edu.duke.cs.osprey.TestBase.*;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;
import java.util.Arrays;

public class TestSimpleGMECFinder {
	
	private static final double EnergyEpsilon = 1e-6;
	
	private static class Problem {
		
		public final SimpleConfSpace confSpace;
		public final ForcefieldParams ffparams;
		public final EnergyCalculator ecalc;
		public final ConfEnergyCalculator confEcalc;
		public final EnergyMatrix emat;
		
		public Problem(Strand ... strands) {
			confSpace = new SimpleConfSpace.Builder().addStrands(strands).build();
			ffparams = new ForcefieldParams();
			ecalc = new EnergyCalculator.Builder(confSpace, ffparams).build();
			confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();
			emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();
		}

		public SimpleGMECFinder makeFinder() {
			return new SimpleGMECFinder.Builder(
				new ConfAStarTree.Builder(emat, confSpace).build(),
				confEcalc
			)
			.build();
		}

		public SimpleGMECFinder makeExternalFinder() {
			return new SimpleGMECFinder.Builder(
				new ConfAStarTree.Builder(emat, confSpace)
					.useExternalMemory()
					.build(),
				confEcalc
			)
			.useExternalMemory()
			.build();
		}

		public SimpleGMECFinder makeConfDBFinder(File confdbFile, Integer interruptAtConfNum) {
			return new SimpleGMECFinder.Builder(
				new ConfAStarTree.Builder(emat, confSpace).build(),
				new ConfEnergyCalculator(confEcalc) {

					int numConfs = 0;

					// override the energy calculation to throw an error after a few confs
					@Override
					public EnergyCalculator.EnergiedParametricMolecule calcEnergy(RCTuple frag, ResidueInteractions inters) {

						if (interruptAtConfNum != null && numConfs++ == interruptAtConfNum) {
							throw new Error("Interrupted!");
						}

						return confEcalc.calcEnergy(frag, inters);
					}
				}
			)
			.setConfDB(confdbFile)
			.build();
		}
	}
	
	private static Problem problemDiscrete;
	private static Problem problemContinuous;
	private static Problem problemMultipleStrands;
	private static Problem problemBigContinuous;
	
	@BeforeClass
	public static void beforeClass() {
		
		Molecule mol = PDBIO.readFile("examples/python.GMEC/1CC8.ss.pdb");
		
		Strand strand1 = new Strand.Builder(mol).build();
		strand1.flexibility.get("A2").setLibraryRotamers("ALA", "GLY");
		strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "VAL");
		strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
		problemDiscrete = new Problem(strand1);
		
		Strand strand2 = new Strand.Builder(mol).setResidues("A2", "A30").build();
		strand2.flexibility.get("A2").setLibraryRotamers("ALA", "GLY");
		strand2.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "VAL", "ARG").setContinuous(10);
		strand2.flexibility.get("A4").addWildTypeRotamers();
		problemContinuous = new Problem(strand2);
		
		Strand strand3 = new Strand.Builder(mol).setResidues("A2", "A30").build();
		strand3.flexibility.get("A2").addWildTypeRotamers();
		strand3.flexibility.get("A3").addWildTypeRotamers();
		strand3.flexibility.get("A4").addWildTypeRotamers();
		Strand strand4 = new Strand.Builder(mol).setResidues("A31", "A60").build();
		strand4.flexibility.get("A31").addWildTypeRotamers();
		strand4.flexibility.get("A32").addWildTypeRotamers();
		strand4.flexibility.get("A33").addWildTypeRotamers();
		problemMultipleStrands = new Problem(strand3, strand4);

		Strand strand5 = new Strand.Builder(mol).setResidues("A2", "A30").build();
		for (String resNum : Arrays.asList("A2", "A3", "A4", "A5", "A6", "A7")) {
			strand5.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType, "ALA", "GLY")
				.setContinuous();
		}
		problemBigContinuous = new Problem(strand5);
	}

	@Test
	public void findDiscrete() {
		EnergiedConf conf = problemDiscrete.makeFinder().find();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-30.705504, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-30.705504, EnergyEpsilon));
	}
	
	@Test
	public void findDiscreteWindowZero() {
		Queue<EnergiedConf> confs = problemDiscrete.makeFinder().find(0);
		assertThat(confs.size(), is(1L));
		
		EnergiedConf conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-30.705504, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-30.705504, EnergyEpsilon));
	}
	
	@Test
	public void findDiscreteWindowOne() {
		Queue<EnergiedConf> confs = problemDiscrete.makeFinder().find(1);
		assertThat(confs.size(), is(4L));
		
		EnergiedConf conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-30.705504, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-30.705504, EnergyEpsilon));
		
		conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 5 }));
		assertThat(conf.getEnergy(), isAbsolutely(-30.241032, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-30.241032, EnergyEpsilon));
		
		conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 6, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-29.981955, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-29.981955, EnergyEpsilon));
		
		conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 7, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-29.748971, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-29.748971, EnergyEpsilon));
	}
	
	@Test
	public void findContinuous() {
		EnergiedConf conf = problemContinuous.makeFinder().find();
		assertThat(conf.getAssignments(), is(new int[] { 1, 26, 0 }));
		assertThat(conf.getEnergy(), isAbsolutely(-38.465807, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-38.566297, EnergyEpsilon));
	}
	
	@Test
	public void findContinuousWindow() {
		Queue<EnergiedConf> confs = problemContinuous.makeFinder().find(0.3);
		assertThat(confs.size(), is(3L));
		
		EnergiedConf conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 26, 0 }));
		assertThat(conf.getEnergy(), isAbsolutely(-38.465807, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-38.566297, EnergyEpsilon));
		
		conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 25, 0 }));
		assertThat(conf.getEnergy(), isAbsolutely(-38.243730, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-38.391590, EnergyEpsilon));
		
		conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 29, 0 }));
		assertThat(conf.getEnergy(), isAbsolutely(-38.166219, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-38.254643, EnergyEpsilon));
	}
	
	@Test
	public void findContinuousWindowExternal() {
		ExternalMemory.use(64, () -> {
			Queue<EnergiedConf> confs = problemContinuous.makeExternalFinder().find(0.3);
			assertThat(confs.size(), is(3L));
			
			EnergiedConf conf = confs.poll();
			assertThat(conf.getAssignments(), is(new int[] { 1, 26, 0 }));
			assertThat(conf.getEnergy(), isAbsolutely(-38.465807, EnergyEpsilon));
			assertThat(conf.getScore(), isAbsolutely(-38.566297, EnergyEpsilon));
			
			conf = confs.poll();
			assertThat(conf.getAssignments(), is(new int[] { 1, 25, 0 }));
			assertThat(conf.getEnergy(), isAbsolutely(-38.243730, EnergyEpsilon));
			assertThat(conf.getScore(), isAbsolutely(-38.391590, EnergyEpsilon));
			
			conf = confs.poll();
			assertThat(conf.getAssignments(), is(new int[] { 1, 29, 0 }));
			assertThat(conf.getEnergy(), isAbsolutely(-38.166219, EnergyEpsilon));
			assertThat(conf.getScore(), isAbsolutely(-38.254643, EnergyEpsilon));
		});
	}

	@Test
	public void findWithResumeNotInterrupted() {

		try (TempFile confDBFile = new TempFile("findgmec.conf.db")) {

			Queue<EnergiedConf> confs = problemBigContinuous.makeConfDBFinder(confDBFile, null).find(0.3);
			assertThat(confs.size(), is(2L));

			EnergiedConf conf = confs.poll();
			assertThat(conf.getAssignments(), is(new int[] { 1, 3, 8, 20, 4, 4 }));
			assertThat(conf.getEnergy(), isAbsolutely(-69.152448, EnergyEpsilon));
			assertThat(conf.getScore(), isAbsolutely(-69.995731, EnergyEpsilon));

			conf = confs.poll();
			assertThat(conf.getAssignments(), is(new int[] { 1, 3, 8, 20, 4, 5 }));
			assertThat(conf.getEnergy(), isAbsolutely(-68.958418, EnergyEpsilon));
			assertThat(conf.getScore(), isAbsolutely(-69.861970, EnergyEpsilon));

			// the conf db should have all 27 confs
			assertThat(confDBFile.exists(), is(true));
			try (ConfDB confdb = new ConfDB(problemBigContinuous.confSpace, confDBFile)) {
				assertThat(confdb.new ConfTable(SimpleGMECFinder.ConfDBTableName).size(), is(27L));
			}
		}
	}

	@Test
	public void findWithResumeInterrupted() {

		try (TempFile confDBFile = new TempFile("findgmec.conf.db")) {

			try {
				// throw an error before the 10th calculated (not cached) energy
				problemBigContinuous.makeConfDBFinder(confDBFile, 10).find(0.3);

				// should throw an error, fail if not
				fail("Expected error, but didn't find one");

			} catch (Error err) {
				// error expected here, keep going
			}

			// there should be a confdb with confs in it
			assertThat(confDBFile.exists(), is(true));
			try (ConfDB confdb = new ConfDB(problemBigContinuous.confSpace, confDBFile)) {
				assertThat(confdb.new ConfTable(SimpleGMECFinder.ConfDBTableName).size(), is(10L));
			}

			// resume the computation using the db
			Queue<EnergiedConf> confs = problemBigContinuous.makeConfDBFinder(confDBFile, null).find(0.3);
			assertThat(confs.size(), is(2L));

			EnergiedConf conf = confs.poll();
			assertThat(conf.getAssignments(), is(new int[] { 1, 3, 8, 20, 4, 4 }));
			assertThat(conf.getEnergy(), isAbsolutely(-69.152448, EnergyEpsilon));
			assertThat(conf.getScore(), isAbsolutely(-69.995731, EnergyEpsilon));

			conf = confs.poll();
			assertThat(conf.getAssignments(), is(new int[] { 1, 3, 8, 20, 4, 5 }));
			assertThat(conf.getEnergy(), isAbsolutely(-68.958418, EnergyEpsilon));
			assertThat(conf.getScore(), isAbsolutely(-69.861970, EnergyEpsilon));

			// the conf db should have all 27 confs
			assertThat(confDBFile.exists(), is(true));
			try (ConfDB confdb = new ConfDB(problemBigContinuous.confSpace, confDBFile)) {
				assertThat(confdb.new ConfTable(SimpleGMECFinder.ConfDBTableName).size(), is(27L));
			}
		}
	}

	@Test
	public void findWithResumeInterruptedTwice() {

		try (TempFile confDBFile = new TempFile("findgmec.conf.db")) {

			try {
				// throw an error before the 10th calculated (not cached) energy
				problemBigContinuous.makeConfDBFinder(confDBFile, 10).find(0.3);

				// should throw an error, fail if not
				fail("Expected error, but didn't find one");

			} catch (Error err) {
				// error expected here, keep going
			}

			// there should be a confdb with confs in it
			assertThat(confDBFile.exists(), is(true));
			try (ConfDB confdb = new ConfDB(problemBigContinuous.confSpace, confDBFile)) {
				assertThat(confdb.new ConfTable(SimpleGMECFinder.ConfDBTableName).size(), is(10L));
			}

			// do another batch
			try {
				// throw an error before the 10th calculated (not cached) energy
				problemBigContinuous.makeConfDBFinder(confDBFile, 10).find(0.3);

				// should throw an error, fail if not
				fail("Expected error, but didn't find one");

			} catch (Error err) {
				// error expected here, keep going
			}

			// there should be a confdb with confs in it
			assertThat(confDBFile.exists(), is(true));
			try (ConfDB confdb = new ConfDB(problemBigContinuous.confSpace, confDBFile)) {
				assertThat(confdb.new ConfTable(SimpleGMECFinder.ConfDBTableName).size(), is(20L));
			}

			// resume the computation using the db
			Queue<EnergiedConf> confs = problemBigContinuous.makeConfDBFinder(confDBFile, null).find(0.3);
			assertThat(confs.size(), is(2L));

			EnergiedConf conf = confs.poll();
			assertThat(conf.getAssignments(), is(new int[] { 1, 3, 8, 20, 4, 4 }));
			assertThat(conf.getEnergy(), isAbsolutely(-69.152448, EnergyEpsilon));
			assertThat(conf.getScore(), isAbsolutely(-69.995731, EnergyEpsilon));

			conf = confs.poll();
			assertThat(conf.getAssignments(), is(new int[] { 1, 3, 8, 20, 4, 5 }));
			assertThat(conf.getEnergy(), isAbsolutely(-68.958418, EnergyEpsilon));
			assertThat(conf.getScore(), isAbsolutely(-69.861970, EnergyEpsilon));

			// the conf db should have all 27 confs
			assertThat(confDBFile.exists(), is(true));
			try (ConfDB confdb = new ConfDB(problemBigContinuous.confSpace, confDBFile)) {
				assertThat(confdb.new ConfTable(SimpleGMECFinder.ConfDBTableName).size(), is(27L));
			}
		}
	}

	@Test
	public void findMultipleStrands() {
		EnergiedConf conf = problemMultipleStrands.makeFinder().find();
		assertThat(conf.getAssignments(), is(new int[] { 0, 0, 0, 0, 0, 0 }));
		assertThat(conf.getEnergy(), isAbsolutely(-84.555275, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-84.555275, EnergyEpsilon));
	}
}
