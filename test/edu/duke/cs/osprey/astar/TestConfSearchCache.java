package edu.duke.cs.osprey.astar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfSearchCache;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;


public class TestConfSearchCache {

	@Test
	public void treeReinstantiation() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A2", "A3", "A4")) {
			strand.flexibility.get(resNum).setLibraryRotamers("VAL");
		}

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		RCs rcs = new RCs(confSpace);

		// calc an emat
		EnergyMatrix emat;
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(8))
			.build()
		) {
			emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();
		}

		// get the full list of expected conformations
		List<ConfSearch.ScoredConf> expectedConfs = new ConfAStarTree.Builder(emat, rcs)
			.setTraditional()
			.build()
			.nextConfs(Double.POSITIVE_INFINITY);
		assertThat(expectedConfs.size(), is(27));

		// make a cache with one tree
			 ConfSearchCache<SimpleConfSpace> cache = new ConfSearchCache<>(1);
			 ConfSearch tree = cache.getOrMake(confSpace, () ->
			new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build()
		);

		assertThat(tree.getNumConformations().intValueExact(), is(expectedConfs.size()));

		for (int i=0; i<10; i++) {
			assertThat(tree.nextConf(), is(expectedConfs.get(i)));
		}

		// clear the trees and force re-instantiation
		cache.clearRefs();

		for (int i=10; i<20; i++) {
			assertThat(tree.nextConf(), is(expectedConfs.get(i)));
		}

		// clear the trees and force re-instantiation again, just for fun
		cache.clearRefs();

		for (int i=20; i<27; i++) {
			assertThat(tree.nextConf(), is(expectedConfs.get(i)));
		}

		// make anything after that is all nulls
		for (int i=27; i<30; i++) {
			assertThat(tree.nextConf(), is(nullValue()));
		}
	}
}
