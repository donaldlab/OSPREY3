package edu.duke.cs.osprey.gmec;

import static edu.duke.cs.osprey.TestBase.*;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.minimization.SimpleConfMinimizer;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

public class TestSimpleGMECFinder {
	
	private static final double EnergyEpsilon = 1e-6;
	
	private static SimpleConfSpace confSpaceDiscrete;
	private static EnergyMatrix ematDiscrete;
	
	private static SimpleConfSpace confSpaceContinuous;
	private static EnergyMatrix ematContinuous;
	private static ConfEnergyCalculator confMinimizer;
	
	@BeforeClass
	public static void beforeClass() {
		
		Molecule mol = PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb");
		
		Strand strand = Strand.builder(mol).build();
		strand.flexibility.get(2).setLibraryRotamers("ALA", "GLY");
		strand.flexibility.get(3).setLibraryRotamers(Strand.WildType, "VAL");
		strand.flexibility.get(4).setLibraryRotamers();
		
		confSpaceDiscrete = SimpleConfSpace.build(strand);
		ematDiscrete = SimplerEnergyMatrixCalculator.build(confSpaceDiscrete).calcEnergyMatrix();
		
		strand = Strand.builder(mol).setResidues(2, 30).build();
		strand.flexibility.get(2).setLibraryRotamers("ALA", "GLY");
		strand.flexibility.get(3).setLibraryRotamers(Strand.WildType, "VAL", "ARG").setContinuous(10);
		strand.flexibility.get(4).addWildTypeRotamers();
		
		confSpaceContinuous = SimpleConfSpace.build(strand);
		ematContinuous = SimplerEnergyMatrixCalculator.build(confSpaceContinuous).calcEnergyMatrix();
		confMinimizer = SimpleConfMinimizer.builder(confSpaceContinuous).build();
	}

	@Test
	public void findDiscrete() {
		EnergiedConf conf = SimpleGMECFinder.builder(confSpaceDiscrete, ematDiscrete).build().find();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-30.705504, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-30.705504, EnergyEpsilon));
	}
	
	@Test
	public void findDiscreteWindowZero() {
		List<EnergiedConf> confs = SimpleGMECFinder.builder(confSpaceDiscrete, ematDiscrete).build().find(0);
		assertThat(confs.size(), is(1));
		
		EnergiedConf conf = confs.get(0);
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 4 }));
		assertThat(conf.getEnergy(), isAbsolutely(-30.705504, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-30.705504, EnergyEpsilon));
	}
	
	@Test
	public void findDiscreteWindowOne() {
		List<EnergiedConf> confs = SimpleGMECFinder.builder(confSpaceDiscrete, ematDiscrete).build().find(1);
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
		EnergiedConf conf = SimpleGMECFinder.builder(confSpaceContinuous, ematContinuous)
			.setEnergyCalculator(confMinimizer)
			.build()
			.find();
		assertThat(conf.getAssignments(), is(new int[] { 1, 26, 0 }));
		assertThat(conf.getEnergy(), isAbsolutely(-38.465807, EnergyEpsilon));
		assertThat(conf.getScore(), isAbsolutely(-38.566297, EnergyEpsilon));
	}
	
	@Test
	public void findContinuousWindow() {
		List<EnergiedConf> confs = SimpleGMECFinder.builder(confSpaceContinuous, ematContinuous)
			.setEnergyCalculator(confMinimizer)
			.build()
			.find(0.1);
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
