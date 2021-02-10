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

package edu.duke.cs.osprey.pruning;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
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


public class TestUnprunedConfBounding {

	private static SimpleConfSpace confSpace;
	private static EnergyMatrix emat;

	@BeforeClass
	public static void beforeClass() {

		// get a conf space small enough we can exhaustively enumerate
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A3").setLibraryRotamers(Strand.WildType);
		strand.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
		strand.flexibility.get("A6").setLibraryRotamers(Strand.WildType);
		confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		// calc the energy matrix
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build())
		{
			emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();
		}
	}

	@Test
	public void fullConfSpace() {

		// get a NOP pruning matrix
		PruningMatrix pmat = new PruningMatrix(confSpace);

		long size = enumerate(pmat);
		assertThat(upper(pmat), greaterThanOrEqualTo(size));
		assertThat(lower(pmat), lessThanOrEqualTo(size));
	}

	@Test
	public void stericPruning() {

		// use only steric pruning (prunes 8 singles, but no pairs)
		PruningMatrix pmat = new SimpleDEE.Runner()
			.setThreshold(100.0)
			.setGoldsteinDiffThreshold(null)
			.run(confSpace, emat);

		long size = enumerate(pmat);
		assertThat(upper(pmat), greaterThanOrEqualTo(size));
		assertThat(lower(pmat), lessThanOrEqualTo(size));
	}

	@Test
	public void intervalPruning5() {

		// use only interval pruning (prunes 14 singles, 2 pairs)
		PruningMatrix pmat = new SimpleDEE.Runner()
			.setThreshold(null)
			.setGoldsteinDiffThreshold(5.0)
			.run(confSpace, emat);

		long size = enumerate(pmat);
		assertThat(upper(pmat), greaterThanOrEqualTo(size));
		assertThat(lower(pmat), lessThanOrEqualTo(size));
	}

	@Test
	public void intervalPruning1() {

		// use only interval pruning (prunes 16 singles, 3 pairs)
		PruningMatrix pmat = new SimpleDEE.Runner()
			.setThreshold(null)
			.setGoldsteinDiffThreshold(1.0)
			.run(confSpace, emat);

		long size = enumerate(pmat);
		assertThat(upper(pmat), greaterThanOrEqualTo(size));
		assertThat(lower(pmat), lessThanOrEqualTo(size));
	}

	@Test
	public void fullCoverPairs() {

		PruningMatrix pmat = new PruningMatrix(confSpace);

		// pos0 -> 8 RCs
		// pos1 -> 7 RCs
		// pos2 -> 8 RCs

		pmat.prunePair(0, 4, 1, 5);
		pmat.prunePair(0, 7, 1, 0);

		pmat.prunePair(0, 3, 2, 4);
		pmat.prunePair(0, 7, 2, 1);

		pmat.prunePair(1, 5, 2, 6);
		pmat.prunePair(1, 4, 2, 0);

		long size = enumerate(pmat);
		assertThat(upper(pmat), greaterThanOrEqualTo(size));
		assertThat(lower(pmat), lessThanOrEqualTo(size));
	}

	@Test
	public void diagonalPairs() {

		PruningMatrix pmat = new PruningMatrix(confSpace);

		pmat.prunePair(0, 4, 1, 5);
		pmat.prunePair(0, 7, 1, 0);

		pmat.prunePair(1, 5, 2, 6);
		pmat.prunePair(1, 4, 2, 0);

		long size = enumerate(pmat);
		assertThat(upper(pmat), greaterThanOrEqualTo(size));
		assertThat(lower(pmat), lessThanOrEqualTo(size));
	}

	@Test
	public void offDiagonalPairs() {

		PruningMatrix pmat = new PruningMatrix(confSpace);

		// pick only "off-diagonal" pairs (e.g., 02)
		pmat.prunePair(0, 3, 2, 4);
		pmat.prunePair(0, 7, 2, 1);

		long size = enumerate(pmat);
		assertThat(upper(pmat), greaterThanOrEqualTo(size));
		assertThat(lower(pmat), lessThanOrEqualTo(size));
	}

	@Test
	public void biggerConfSpace() {

		// get a conf space small enough we can exhaustively enumerate
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A3", "A4", "A5", "A6")) {
			strand.flexibility.get(resNum).setLibraryRotamers(Strand.WildType, "ARG", "LYS");
		}
		confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		// calc the energy matrix
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build())
		{
			emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();

			// use only interval pruning
			PruningMatrix pmat = new SimpleDEE.Runner()
				.setThreshold(null)
				.setGoldsteinDiffThreshold(10.0)
				.run(confSpace, emat);

			long size = enumerate(pmat);
			assertThat(upper(pmat), greaterThanOrEqualTo(size));
			assertThat(lower(pmat), lessThanOrEqualTo(size));
		}
	}

	@Test
	public void biggerConfSpace2() {

		// get a conf space small enough we can exhaustively enumerate
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A3", "A4", "A5", "A6")) {
			strand.flexibility.get(resNum).setLibraryRotamers(Strand.WildType, "VAL", "LEU");
		}
		confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		// calc the energy matrix
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build())
		{
			emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();

			// use only interval pruning
			PruningMatrix pmat = new SimpleDEE.Runner()
				.setThreshold(null)
				.setGoldsteinDiffThreshold(10.0)
				.run(confSpace, emat);

			long size = enumerate(pmat);
			assertThat(upper(pmat), greaterThanOrEqualTo(size));
			assertThat(lower(pmat), lessThanOrEqualTo(size));
		}
	}

	private static long enumerate(PruningMatrix pmat) {
		ConfAStarTree astar = new ConfAStarTree.Builder(emat, pmat)
			.setTraditional()
			.build();
		long size = 0L;
		while (astar.nextConf() != null) {
			size++;
		}
		return size;
	}

	private static long upper(PruningMatrix pmat) {
		return pmat.calcUnprunedConfsUpperBound().longValueExact();
	}

	private static long lower(PruningMatrix pmat) {
		return pmat.calcUnprunedConfsLowerBound().longValueExact();
	}
}
