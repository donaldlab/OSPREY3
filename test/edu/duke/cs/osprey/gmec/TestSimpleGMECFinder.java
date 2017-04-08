package edu.duke.cs.osprey.gmec;

import static edu.duke.cs.osprey.TestBase.*;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.MinimizingConfEnergyCalculator;
import edu.duke.cs.osprey.energy.MinimizingFragmentEnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

public class TestSimpleGMECFinder {
	
	private static final double EnergyEpsilon = 1e-6;
	
	private static class Problem {
		
		public final SimpleConfSpace confSpace;
		public final ForcefieldParams ffparams;
		public final MinimizingFragmentEnergyCalculator ecalc;
		public final EnergyMatrix emat;
		
		public Problem(Strand strand) {
			confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
			ffparams = new ForcefieldParams();
			ecalc = new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams).build();
			emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();
		}
		
		public SimpleGMECFinder makeGMECFinder() {
			return new SimpleGMECFinder.Builder(
				confSpace,
				new ConfAStarTree.Builder(emat, confSpace).build(),
				new MinimizingConfEnergyCalculator.Builder(ecalc).build()
			).build();
		}
	}
	
	private static Problem problemDiscrete;
	private static Problem problemContinuous;
	
	@BeforeClass
	public static void beforeClass() {
		
		Molecule mol = PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb");
		
		Strand strand = new Strand.Builder(mol).build();
		strand.flexibility.get(2).setLibraryRotamers("ALA", "GLY");
		strand.flexibility.get(3).setLibraryRotamers(Strand.WildType, "VAL");
		strand.flexibility.get(4).setLibraryRotamers();
		problemDiscrete = new Problem(strand);
		
		strand = new Strand.Builder(mol).setResidues(2, 30).build();
		strand.flexibility.get(2).setLibraryRotamers("ALA", "GLY");
		strand.flexibility.get(3).setLibraryRotamers(Strand.WildType, "VAL", "ARG").setContinuous(10);
		strand.flexibility.get(4).addWildTypeRotamers();
		problemContinuous = new Problem(strand);
	}

	@Test
	public void findDiscrete() {
		EnergiedConf conf = problemDiscrete.makeGMECFinder().find();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-30.705504, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-30.705504, EnergyEpsilon));
	}
	
	@Test
	public void findDiscreteWindowZero() {
		List<EnergiedConf> confs = problemDiscrete.makeGMECFinder().find(0);
		assertThat(confs.size(), is(1));
		
		EnergiedConf conf = confs.get(0);
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-30.705504, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-30.705504, EnergyEpsilon));
	}
	
	@Test
	public void findDiscreteWindowOne() {
		List<EnergiedConf> confs = problemDiscrete.makeGMECFinder().find(1);
		assertThat(confs.size(), is(4));
		
		EnergiedConf conf = confs.get(0);
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-30.705504, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-30.705504, EnergyEpsilon));
		
		conf = confs.get(1);
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 5 }));
		assertThat(conf.getEnergy(), isAbsolutely(-30.241032, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-30.241032, EnergyEpsilon));
		
		conf = confs.get(2);
		assertThat(conf.getAssignments(), is(new int[] { 1, 6, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-29.981955, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-29.981955, EnergyEpsilon));
		
		conf = confs.get(3);
		assertThat(conf.getAssignments(), is(new int[] { 1, 7, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-29.748971, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-29.748971, EnergyEpsilon));
	}
	
	@Test
	public void findContinuous() {
		EnergiedConf conf = problemContinuous.makeGMECFinder().find();
		assertThat(conf.getAssignments(), is(new int[] { 1, 26, 0 }));
		assertThat(conf.getEnergy(), isAbsolutely(-38.465807, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-38.566297, EnergyEpsilon));
	}
	
	@Test
	public void findContinuousWindow() {
		List<EnergiedConf> confs = problemContinuous.makeGMECFinder().find(0.1);
		assertThat(confs.size(), is(3));
		
		EnergiedConf conf = confs.get(0);
		assertThat(conf.getAssignments(), is(new int[] { 1, 26, 0 }));
		assertThat(conf.getEnergy(), isAbsolutely(-38.465807, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-38.566297, EnergyEpsilon));
		
		conf = confs.get(1);
		assertThat(conf.getAssignments(), is(new int[] { 1, 25, 0 }));
		assertThat(conf.getEnergy(), isAbsolutely(-38.243730, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-38.391590, EnergyEpsilon));
		
		conf = confs.get(2);
		assertThat(conf.getAssignments(), is(new int[] { 1, 21, 0 }));
		assertThat(conf.getEnergy(), isAbsolutely(-37.755997, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-38.695671, EnergyEpsilon));
	}
}
