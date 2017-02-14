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
import edu.duke.cs.osprey.structure.PDBIO;

public class TestSimpleGMECFinder {
	
	private static final double EnergyEpsilon = 1e-6;
	
	private static SimpleConfSpace confSpace;
	private static EnergyMatrix emat;
	
	@BeforeClass
	public static void beforeClass() {
		
		Strand strand = Strand.builder(PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb")).build();
		strand.flexibility.get(2).setLibraryRotamers("ALA", "GLY");
		strand.flexibility.get(3).setLibraryRotamers(Strand.WildType, "VAL");
		strand.flexibility.get(4).setLibraryRotamers();
		
		confSpace = SimpleConfSpace.build(strand);
		emat = SimplerEnergyMatrixCalculator.build(confSpace).calcEnergyMatrix();
	}

	@Test
	public void find() {
		EnergiedConf conf = SimpleGMECFinder.builder(confSpace, emat).build().find();
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 4 }));
		assertThat(conf.getScore(), isAbsolutely(-30.705504, EnergyEpsilon));
		assertThat(conf.getEnergy(), isAbsolutely(-30.705504, EnergyEpsilon));
	}
	
	@Test
	public void findWindowZero() {
		List<EnergiedConf> confs = SimpleGMECFinder.builder(confSpace, emat).build().find(0);
		assertThat(confs.size(), is(1));
		
		EnergiedConf conf = confs.get(0);
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 4 }));
		assertThat(conf.getScore(), isAbsolutely(-30.705504, EnergyEpsilon));
		assertThat(conf.getEnergy(), isAbsolutely(-30.705504, EnergyEpsilon));
	}
	
	@Test
	public void findWindowOne() {
		List<EnergiedConf> confs = SimpleGMECFinder.builder(confSpace, emat).build().find(1);
		assertThat(confs.size(), is(4));
		
		EnergiedConf conf = confs.get(0);
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 4 }));
		assertThat(conf.getScore(), isAbsolutely(-30.705504, EnergyEpsilon));
		assertThat(conf.getEnergy(), isAbsolutely(-30.705504, EnergyEpsilon));
		
		conf = confs.get(1);
		assertThat(conf.getAssignments(), is(new int[] { 1, 3, 5 }));
		assertThat(conf.getScore(), isAbsolutely(-30.241032, EnergyEpsilon));
		assertThat(conf.getEnergy(), isAbsolutely(-30.241032, EnergyEpsilon));
		
		conf = confs.get(2);
		assertThat(conf.getAssignments(), is(new int[] { 1, 6, 4 }));
		assertThat(conf.getScore(), isAbsolutely(-29.981955, EnergyEpsilon));
		assertThat(conf.getEnergy(), isAbsolutely(-29.981955, EnergyEpsilon));
		
		conf = confs.get(3);
		assertThat(conf.getAssignments(), is(new int[] { 1, 7, 4 }));
		assertThat(conf.getScore(), isAbsolutely(-29.748971, EnergyEpsilon));
		assertThat(conf.getEnergy(), isAbsolutely(-29.748971, EnergyEpsilon));
	}
}
