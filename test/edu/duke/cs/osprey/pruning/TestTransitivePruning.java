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

import static edu.duke.cs.osprey.tools.Log.log;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.Test;

import java.util.*;


public class TestTransitivePruning {

	@Test
	public void test1CC8_1() {

		Strand design = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A16")) {
			design.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		}

		check(design);
	}

	@Test
	public void test1CC8_2() {

		Strand design = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A16", "A17")) {
			design.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		}

		check(design);
	}

	@Test
	public void test1CC8_3() {

		Strand design = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A16", "A17", "A18")) {
			design.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		}

		check(design);
	}

	@Test
	public void test1CC8_4() {

		Strand design = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A16", "A17", "A18", "A19")) {
			design.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		}

		check(design);
	}

	@Test
	public void test1CC8_5() {

		Strand design = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A16", "A17", "A18", "A19", "A20")) {
			design.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		}

		check(design);
	}

	@Test
	public void test1CC8_6() {

		Strand design = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A16", "A17", "A18", "A19", "A20", "A21")) {
			design.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		}

		check(design);
	}

	@Test
	public void test1CC8_7() {

		Strand design = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A16", "A17", "A18", "A19", "A20", "A21", "A22")) {
			design.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		}

		check(design);
	}

	@Test
	public void test1CC8_8() {

		Strand design = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A16", "A17", "A18", "A19", "A20", "A21", "A22", "A23")) {
			design.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		}

		check(design);
	}

	@Test
	public void test1CC8_9() {

		Strand design = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A16", "A17", "A18", "A19", "A20", "A21", "A22", "A23", "A24")) {
			design.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		}

		check(design);
	}

	@Test
	public void test1CC8_10() {

		Strand design = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A16", "A17", "A18", "A19", "A20", "A21", "A22", "A23", "A24", "A25")) {
			design.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers();
		}

		check(design);
	}

	private void check(Strand strand) {

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrands(strand).build();

		// calc the pruning matrix
		PruningMatrix pmat;
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.build()
				.calcEnergyMatrix();

			pmat = new SimpleDEE.Runner()
				.setSinglesThreshold(100.0)
				.setPairsThreshold(100.0)
				.setGoldsteinDiffThreshold(5.0)
				.run(confSpace, emat);
		}

		check(confSpace, pmat);
	}

	private void check(SimpleConfSpace confSpace, PruningMatrix pmat) {

		// if a tuple cannot be part of a full conformation due to pruning
		// (ie, all paths are pruned before reaching a leaf node),
		// then we expect that tuple should be pruned by transitive pruning

		// find out which tuples should be transitively pruned by brute force

		Set<RCTuple> prunedSingles = new HashSet<>();
		Set<RCTuple> unprunedSingles = new HashSet<>();
		pmat.forEachUnprunedSingle((pos1, rc1) -> {
			RCTuple tuple = new RCTuple(pos1, rc1).sorted();
			if (hasLeaf(confSpace, pmat, tuple)) {
				unprunedSingles.add(tuple);
			} else {
				prunedSingles.add(tuple);
			}
			return PruningMatrix.IteratorCommand.Continue;
		});


		Set<RCTuple> prunedPairs = new HashSet<>();
		Set<RCTuple> unprunedPairs = new HashSet<>();
		pmat.forEachUnprunedPair((pos1, rc1, pos2, rc2) -> {
			RCTuple tuple = new RCTuple(pos1, rc1, pos2, rc2).sorted();
			if (hasLeaf(confSpace, pmat, tuple)) {
				unprunedPairs.add(tuple);
			} else {
				prunedPairs.add(tuple);
			}
			return PruningMatrix.IteratorCommand.Continue;
		});

		Set<RCTuple> prunedTriples = new HashSet<>();
		Set<RCTuple> unprunedTriples = new HashSet<>();
		pmat.forEachUnprunedTriple((pos1, rc1, pos2, rc2, pos3, rc3) -> {
			RCTuple tuple = new RCTuple(pos1, rc1, pos2, rc2, pos3, rc3).sorted();
			if (hasLeaf(confSpace, pmat, tuple)) {
				unprunedTriples.add(tuple);
			} else {
				prunedTriples.add(tuple);
			}
			return PruningMatrix.IteratorCommand.Continue;
		});

		// TEMP
		log("confs: %s", new RCs(confSpace).getNumConformations());
		log("singles: %d,%d   pairs: %d,%d   triples: %d,%d",
			prunedSingles.size(), unprunedSingles.size(),
			prunedPairs.size(), unprunedPairs.size(),
			prunedTriples.size(), unprunedTriples.size()
		);

		// all the expected tuples shouldn't be pruned yet
		for (RCTuple tuple : prunedSingles) {
			assertThat(pmat.isPruned(tuple), is(false));
		}
		for (RCTuple tuple : prunedPairs) {
			assertThat(pmat.isPruned(tuple), is(false));
		}
		for (RCTuple tuple : prunedTriples) {
			assertThat(pmat.isPruned(tuple), is(false));
		}

		// use the transitive pruner
		TransitivePruner pruner = new TransitivePruner(confSpace);
		TaskExecutor tasks = new TaskExecutor();
		pruner.pruneSingles(pmat, tasks);
		pruner.prunePairs(pmat, tasks);
		pruner.pruneTriples(pmat, tasks);

		// all the expected tuples should be pruned
		for (RCTuple tuple : prunedSingles) {
			assertThat(pmat.isPruned(tuple), is(true));
		}
		for (RCTuple tuple : unprunedSingles) {
			assertThat(pmat.isPruned(tuple), is(false));
		}
		for (RCTuple tuple : prunedPairs) {
			assertThat(pmat.isPruned(tuple), is(true));
		}
		for (RCTuple tuple : unprunedPairs) {
			assertThat(pmat.isPruned(tuple), is(false));
		}
		for (RCTuple tuple : prunedTriples) {
			assertThat(pmat.isPruned(tuple), is(true));
		}
		for (RCTuple tuple : unprunedTriples) {
			assertThat(pmat.isPruned(tuple), is(false));
		}
	}

	private boolean hasLeaf(SimpleConfSpace confSpace, PruningMatrix pmat, RCTuple tuple) {

		assert (!pmat.isPruned(tuple));

		// collect all the positions: assigned first, then unassigned
		List<SimpleConfSpace.Position> positions = new ArrayList<>();
		for (int pos : tuple.pos) {
			positions.add(confSpace.positions.get(pos));
		}
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			if (!positions.contains(pos)) {
				positions.add(pos);
			}
		}

		// sort unassigned positions by number of unpruned RCs,
		// so we're likely to dead-end faster
		positions.subList(tuple.size(), positions.size())
			.sort(Comparator.comparing(pos -> pos.resConfs.stream()
				.filter(rc -> !pmat.isSinglePruned(pos.index, rc.index))
				.count()
			));

		// do DFS
		Stack<ConfIndex> stack = new Stack<>();

		// start with the DFS node representing the tuple
		ConfIndex root = new ConfIndex(positions.size());
		root.numDefined = tuple.size();
		for (int i=0; i<tuple.size(); i++) {
			root.definedPos[i] = tuple.pos.get(i);
			root.definedRCs[i] = tuple.RCs.get(i);
		}
		root.sortDefined();
		root.updateUndefined();
		stack.push(root);

		while (!stack.isEmpty()) {

			ConfIndex node = stack.pop();

			// hit a leaf node? we're done here
			if (node.numDefined == positions.size()) {
				return true;
			}

			// otherwise, expand the next pos and RCs
			SimpleConfSpace.Position pos = positions.get(node.numDefined);
			for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {

				// if this child was pruned by the pruning matrix, then skip it
				if (isPruned(pmat, node, pos.index, rc.index)) {
					continue;
				}

				// otherwise, expand it
				stack.push(node.assign(pos.index, rc.index));
			}
		}

		return false;
	}

	private boolean isPruned(PruningMatrix pmat, ConfIndex confIndex, int nextPos, int nextRc) {

		// check single
		if (pmat.getOneBody(nextPos, nextRc)) {
			return true;
		}

		// check pairs
		for (int i=0; i<confIndex.numDefined; i++) {
			int pos = confIndex.definedPos[i];
			int rc = confIndex.definedRCs[i];
			assert (pos != nextPos || rc != nextRc);
			if (pmat.getPairwise(pos, rc, nextPos, nextRc)) {
				return true;
			}
		}

		// check triples
		if (pmat.hasHigherOrderTuples()) {

			RCTuple tuple = new RCTuple(0, 0, 0, 0, 0, 0);

			for (int i1=0; i1<confIndex.numDefined; i1++) {
				int pos1 = confIndex.definedPos[i1];
				int rc1 = confIndex.definedRCs[i1];
				assert (pos1 != nextPos || rc1 != nextRc);

				for (int i2=0; i2<i1; i2++) {
					int pos2 = confIndex.definedPos[i2];
					int rc2 = confIndex.definedRCs[i2];
					assert (pos2 != nextPos || rc2 != nextRc);

					tuple.set(pos1, rc1, pos2, rc2, nextPos, nextRc);
					tuple.sortPositions();

					if (pmat.getTuple(tuple)) {
						return true;
					}
				}
			}
		}

		return false;
	}
}
