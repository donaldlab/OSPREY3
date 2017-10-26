package edu.duke.cs.osprey.gmec;

//Making sure DEEGMECFinder results match SimpleGMECFinder


import static edu.duke.cs.osprey.TestBase.*;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

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

public class TestDEEGMECFinder {
	
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
		
		public DEEGMECFinder.Builder makeBuilder() {
			return new DEEGMECFinder.Builder(
                                emat, confSpace, ecalc, confEcalc, null
			);
		}
	}
	
	private static Problem problemDiscrete;
	private static Problem problemContinuous;
	private static Problem problemMultipleStrands;
	
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
	}

	@Test
	public void findDiscrete() {
		EnergiedConf conf = problemDiscrete.makeBuilder().build().calcGMEC();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-30.705504, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-30.705504, EnergyEpsilon));
	}
	
	@Test
	public void findDiscreteWindowZero() {
		Queue<EnergiedConf> confs = problemDiscrete.makeBuilder().build().calcGMEC(0);
		assertThat(confs.size(), is(1L));
		
		EnergiedConf conf = confs.poll();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-30.705504, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-30.705504, EnergyEpsilon));
	}
	
	@Test
	public void findDiscreteWindowOne() {
                DEEGMECFinder gf = problemDiscrete.makeBuilder().build();
                gf.Ew = 1;
		Queue<EnergiedConf> confs = gf.calcGMEC(0);
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
		EnergiedConf conf = problemContinuous.makeBuilder().build().calcGMEC();
		assertThat(conf.getAssignments(), is(new int[] { 1, 26, 0 }));
		assertThat(conf.getEnergy(), isAbsolutely(-38.465807, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-38.566297, EnergyEpsilon));
	}
	
	@Test
	public void findContinuousWindow() {
                DEEGMECFinder gf = problemContinuous.makeBuilder().build();
                gf.Ew = 0.3;
		Queue<EnergiedConf> confs = gf.calcGMEC(0);
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
                        DEEGMECFinder gf = ((DEEGMECFinder.Builder)problemContinuous.makeBuilder()
				.useExternalMemory())
				.build();
                        gf.Ew = 0.3;
			Queue<EnergiedConf> confs = gf.calcGMEC(0);
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
	public void findMultipleStrands() {
		EnergiedConf conf = problemMultipleStrands.makeBuilder().build().calcGMEC();
		assertThat(conf.getAssignments(), is(new int[] { 0, 0, 0, 0, 0, 0 }));
		assertThat(conf.getEnergy(), isAbsolutely(-84.555275, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-84.555275, EnergyEpsilon));
	}
}
