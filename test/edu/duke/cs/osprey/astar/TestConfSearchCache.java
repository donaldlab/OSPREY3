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
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;


public class TestConfSearchCache {

	private static RCs rcs;
	private static EnergyMatrix emat;


	@BeforeClass
	public static void beforeClass() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A2", "A3", "A4")) {
			strand.flexibility.get(resNum).setLibraryRotamers("VAL");
		}

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		rcs = new RCs(confSpace);

		// calc an emat
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(8))
			.build()
		) {
			emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();
		}
	}

	@Test
	public void treeReinstantiation() {

		// get the full list of expected conformations
		List<ConfSearch.ScoredConf> expectedConfs = new ConfAStarTree.Builder(emat, rcs)
			.setTraditional()
			.build()
			.nextConfs(Double.POSITIVE_INFINITY);
		assertThat(expectedConfs.size(), is(27));

		// make a cache with one tree
		ConfSearchCache cache = new ConfSearchCache(1);
		ConfSearchCache.Entry tree = cache.make(() ->
			new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build()
		);

		assertThat(tree.getNumConformations().intValueExact(), is(expectedConfs.size()));

		for (int i=0; i<10; i++) {
			assertThat(tree.nextConf(), is(expectedConfs.get(i)));
		}

		// clear the trees and force re-instantiation
		tree.clearRefs();

		for (int i=10; i<20; i++) {
			assertThat(tree.nextConf(), is(expectedConfs.get(i)));
		}

		// clear the trees and force re-instantiation again, just for fun
		tree.clearRefs();

		for (int i=20; i<27; i++) {
			assertThat(tree.nextConf(), is(expectedConfs.get(i)));
		}

		// make anything after that is all nulls
		for (int i=27; i<30; i++) {
			assertThat(tree.nextConf(), is(nullValue()));
		}
	}

	@Test
	public void unrestrictedCapacity() {

		ConfSearchCache cache = new ConfSearchCache(null);

		ConfSearchCache.Entry tree1 = cache.make(() ->
			new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build()
		);
		assertThat(tree1.isProtected(), is(true));

		ConfSearchCache.Entry tree2 = cache.make(() ->
			new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build()
		);
		assertThat(tree1.isProtected(), is(true));
		assertThat(tree2.isProtected(), is(true));

		ConfSearchCache.Entry tree3 = cache.make(() ->
			new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build()
		);
		assertThat(tree1.isProtected(), is(true));
		assertThat(tree2.isProtected(), is(true));
		assertThat(tree3.isProtected(), is(true));
	}

	@Test
	public void restrictedCapacity() {

		ConfSearchCache cache = new ConfSearchCache(2);

		ConfSearchCache.Entry tree1 = cache.make(() ->
			new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build()
		);
		assertThat(tree1.isProtected(), is(true));

		ConfSearchCache.Entry tree2 = cache.make(() ->
			new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build()
		);
		assertThat(tree1.isProtected(), is(true));
		assertThat(tree2.isProtected(), is(true));

		ConfSearchCache.Entry tree3 = cache.make(() ->
			new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build()
		);
		assertThat(tree1.isProtected(), is(false));
		assertThat(tree2.isProtected(), is(true));
		assertThat(tree3.isProtected(), is(true));

		tree1.nextConf();
		assertThat(tree1.isProtected(), is(true));
		assertThat(tree2.isProtected(), is(false));
		assertThat(tree3.isProtected(), is(true));

		tree2.nextConf();
		assertThat(tree1.isProtected(), is(true));
		assertThat(tree2.isProtected(), is(true));
		assertThat(tree3.isProtected(), is(false));

		tree1.nextConf();
		assertThat(tree1.isProtected(), is(true));
		assertThat(tree2.isProtected(), is(true));
		assertThat(tree3.isProtected(), is(false));

		tree2.nextConf();
		assertThat(tree1.isProtected(), is(true));
		assertThat(tree2.isProtected(), is(true));
		assertThat(tree3.isProtected(), is(false));

		tree1.nextConf();
		assertThat(tree1.isProtected(), is(true));
		assertThat(tree2.isProtected(), is(true));
		assertThat(tree3.isProtected(), is(false));

		tree3.nextConf();
		assertThat(tree1.isProtected(), is(true));
		assertThat(tree2.isProtected(), is(false));
		assertThat(tree3.isProtected(), is(true));
	}
}
