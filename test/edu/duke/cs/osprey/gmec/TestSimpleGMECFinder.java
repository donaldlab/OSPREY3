package edu.duke.cs.osprey.gmec;

import static edu.duke.cs.osprey.TestBase.*;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

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

		public SimpleGMECFinder makeResumeFinder(File logFile, Integer interruptAtConfNum) {
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
			.setResumeLog(logFile)
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

		File logFile = new File("resume.log");
		logFile.delete();
		assertThat(logFile.exists(), is(false));

		Queue<EnergiedConf> confs = problemBigContinuous.makeResumeFinder(logFile, null).find(0.3);
		assertThat(confs.size(), is(2L));

		EnergiedConf conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 8, 20, 4, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-69.152448, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-69.995731, EnergyEpsilon));

		conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 8, 20, 4, 5 }));
		assertThat(conf.getEnergy(), isAbsolutely(-68.958418, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-69.861970, EnergyEpsilon));

		// there should be no resume log
		assertThat(logFile.exists(), is(false));
	}

	@Test
	public void findWithResumeInterrupted() {

		// delete any previous resume log, if any, just in case
		File logFile = new File("resume.log");
		logFile.delete();
		assertThat(logFile.exists(), is(false));

		try {
			// NOTE: the resume logger batches writes, so calc at least 8 conformations
			problemBigContinuous.makeResumeFinder(logFile, 10).find(0.3);

			// should throw an error, fail if not
			fail("Expected error, but didn't find one");

		} catch (Error err) {
			// error expected here, keep going
		}

		// there should be a resume log
		assertThat(logFile.exists(), is(true));

		// resume the computation using the log
		Queue<EnergiedConf> confs = problemBigContinuous.makeResumeFinder(logFile, null).find(0.3);
		assertThat(confs.size(), is(2L));

		EnergiedConf conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 8, 20, 4, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-69.152448, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-69.995731, EnergyEpsilon));

		conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 8, 20, 4, 5 }));
		assertThat(conf.getEnergy(), isAbsolutely(-68.958418, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-69.861970, EnergyEpsilon));

		// the resume log should have been cleaned up
		assertThat(logFile.exists(), is(false));
	}

	@Test
	public void findWithResumeInterruptedTwice() {

		// delete any previous resume log, if any, just in case
		File logFile = new File("resume.log");
		logFile.delete();
		assertThat(logFile.exists(), is(false));

		try {
			// NOTE: the resume logger batches writes, so calc at least 8 conformations
			problemBigContinuous.makeResumeFinder(logFile, 10).find(0.3);

			// should throw an error, fail if not
			fail("Expected error, but didn't find one");

		} catch (Error err) {
			// error expected here, keep going
		}

		// there should be a resume log
		assertThat(logFile.exists(), is(true));

		// do another batch
		try {
			problemBigContinuous.makeResumeFinder(logFile, 10).find(0.3);

			// should throw an error, fail if not
			fail("Expected error, but didn't find one");

		} catch (Error err) {
			// error expected here, keep going
		}

		// resume the computation using the log
		Queue<EnergiedConf> confs = problemBigContinuous.makeResumeFinder(logFile, null).find(0.3);
		assertThat(confs.size(), is(2L));

		EnergiedConf conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 8, 20, 4, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-69.152448, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-69.995731, EnergyEpsilon));

		conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 8, 20, 4, 5 }));
		assertThat(conf.getEnergy(), isAbsolutely(-68.958418, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-69.861970, EnergyEpsilon));

		// the resume log should have been cleaned up
		assertThat(logFile.exists(), is(false));
	}

	@Test
	public void findMultipleStrands() {
		EnergiedConf conf = problemMultipleStrands.makeFinder().find();
		assertThat(conf.getAssignments(), is(new int[] { 0, 0, 0, 0, 0, 0 }));
		assertThat(conf.getEnergy(), isAbsolutely(-84.555275, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-84.555275, EnergyEpsilon));
	}
}
