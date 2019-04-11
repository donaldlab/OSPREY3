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

package edu.duke.cs.osprey.sofea;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.TestBase.TempFile;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.BigExp;
import org.hamcrest.BaseMatcher;
import org.hamcrest.Description;
import org.hamcrest.Matcher;
import org.junit.Test;

import java.util.Arrays;


public class TestFringeDB {

	@Test
	public void sizes() {
		MultiStateConfSpace confSpace = makeConfSpace();

		try (TempFile file = new TempFile("fringe.db")) {

			// create a new FringeDB
			try (FringeDB db = FringeDB.create(confSpace, file, 1024)) {
				assertThat(file.length(), is(1024L));
				assertThat(db.getNumNodes(), is(0L));
				assertThat(db.getCapacity(), is(50L));
			}

			// re-open the existing FringeDB
			try (FringeDB db = FringeDB.open(confSpace, file)) {
				assertThat(file.length(), is(1024L));
				assertThat(db.getNumNodes(), is(0L));
				assertThat(db.getCapacity(), is(50L));
			}
		}
	}

	@Test
	public void rootsSweepConsume() {

		MultiStateConfSpace confSpace = makeConfSpace();

		try (TempFile file = new TempFile("fringe.db")) {
			try (FringeDB db = FringeDB.create(confSpace, file, 1024)) {

				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(Double.NaN)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(Double.NaN)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(Double.NaN)));

				// add the root nodes
				FringeDB.Transaction tx = db.transaction();
				tx.writeRootNode(confSpace.states.get(0), new BigExp(1024.5));
				tx.writeRootNode(confSpace.states.get(1), new BigExp(10.4));
				tx.writeRootNode(confSpace.states.get(2), new BigExp(7.3));

				assertThat(tx.dbHasRoomForCommit(), is(true));
				tx.commit();
				db.finishStep();

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(1024.5)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(10.4)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(7.3)));
				assertThat(db.getNumNodes(), is(3L));

				// do a sweep and consume all the nodes
				tx = db.transaction();
				assertThat(tx.numNodesToRead(), is(3L));

				tx.readNode();
				assertThat(tx.numNodesToRead(), is(2L));
				assertThat(tx.state().index, is(0));
				assertThat(tx.conf(), is(conf(-1, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(1024.5)));

				tx.readNode();
				assertThat(tx.numNodesToRead(), is(1L));
				assertThat(tx.state().index, is(1));
				assertThat(tx.conf(), is(conf(-1, -1, -1, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(10.4)));

				tx.readNode();
				assertThat(tx.numNodesToRead(), is(0L));
				assertThat(tx.state().index, is(2));
				assertThat(tx.conf(), is(conf(-1, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(7.3)));

				tx.commit();
				db.finishStep();

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(Double.NaN)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(Double.NaN)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(Double.NaN)));
				assertThat(db.getNumNodes(), is(0L));
			}
		}
	}

	@Test
	public void rootsSweepKeep() {

		MultiStateConfSpace confSpace = makeConfSpace();

		try (TempFile file = new TempFile("fringe.db")) {
			try (FringeDB db = FringeDB.create(confSpace, file, 1024)) {

				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(Double.NaN)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(Double.NaN)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(Double.NaN)));

				// add the root nodes
				FringeDB.Transaction tx = db.transaction();
				tx.writeRootNode(confSpace.states.get(0), new BigExp(1024.5));
				tx.writeRootNode(confSpace.states.get(1), new BigExp(10.4));
				tx.writeRootNode(confSpace.states.get(2), new BigExp(7.3));

				assertThat(tx.dbHasRoomForCommit(), is(true));
				tx.commit();
				db.finishStep();

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(1024.5)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(10.4)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(7.3)));
				assertThat(db.getNumNodes(), is(3L));

				// do a sweep and keep all the nodes
				tx = db.transaction();

				tx.readNode();
				assertThat(tx.state().index, is(0));
				assertThat(tx.conf(), is(conf(-1, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(1024.5)));
				tx.writeReplacementNode(tx.state(), tx.conf(), tx.zSumUpper());

				tx.readNode();
				assertThat(tx.state().index, is(1));
				assertThat(tx.conf(), is(conf(-1, -1, -1, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(10.4)));
				tx.writeReplacementNode(tx.state(), tx.conf(), tx.zSumUpper());

				tx.readNode();
				assertThat(tx.state().index, is(2));
				assertThat(tx.conf(), is(conf(-1, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(7.3)));
				tx.writeReplacementNode(tx.state(), tx.conf(), tx.zSumUpper());

				tx.commit();
				db.finishStep();

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(1024.5)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(10.4)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(7.3)));
				assertThat(db.getNumNodes(), is(3L));
			}
		}
	}

	@Test
	public void rootsResumeSweepConsume() {

		MultiStateConfSpace confSpace = makeConfSpace();

		try (TempFile file = new TempFile("fringe.db")) {

			// make a new DB with root nodes
			try (FringeDB db = FringeDB.create(confSpace, file, 1024)) {
				FringeDB.Transaction tx = db.transaction();
				tx.writeRootNode(confSpace.states.get(0), new BigExp(1024.5));
				tx.writeRootNode(confSpace.states.get(1), new BigExp(10.4));
				tx.writeRootNode(confSpace.states.get(2), new BigExp(7.3));
				tx.commit();
				db.finishStep();
			}

			// open the DB and consume all the nodes
			try (FringeDB db = FringeDB.open(confSpace, file)) {

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(1024.5)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(10.4)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(7.3)));
				assertThat(db.getNumNodes(), is(3L));

				FringeDB.Transaction tx = db.transaction();

				tx.readNode();
				assertThat(tx.state().index, is(0));
				assertThat(tx.conf(), is(conf(-1, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(1024.5)));

				tx.readNode();
				assertThat(tx.state().index, is(1));
				assertThat(tx.conf(), is(conf(-1, -1, -1, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(10.4)));

				tx.readNode();
				assertThat(tx.state().index, is(2));
				assertThat(tx.conf(), is(conf(-1, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(7.3)));

				tx.commit();
				db.finishStep();

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(Double.NaN)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(Double.NaN)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(Double.NaN)));
				assertThat(db.getNumNodes(), is(0L));
			}

			// we should have an empty DB again
			try (FringeDB db = FringeDB.open(confSpace, file)) {
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(Double.NaN)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(Double.NaN)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(Double.NaN)));
				assertThat(db.getNumNodes(), is(0L));
			}
		}
	}

	@Test
	public void rootsSweepResume() {

		MultiStateConfSpace confSpace = makeConfSpace();

		try (TempFile file = new TempFile("fringe.db")) {

			// make a new DB with root nodes
			try (FringeDB db = FringeDB.create(confSpace, file, 1024)) {
				FringeDB.Transaction tx = db.transaction();
				tx.writeRootNode(confSpace.states.get(0), new BigExp(1024.5));
				tx.writeRootNode(confSpace.states.get(1), new BigExp(10.4));
				tx.writeRootNode(confSpace.states.get(2), new BigExp(7.3));
				tx.commit();
				db.finishStep();
			}

			// open the DB and consume one node
			try (FringeDB db = FringeDB.open(confSpace, file)) {

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(1024.5)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(10.4)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(7.3)));
				assertThat(db.getNumNodes(), is(3L));

				FringeDB.Transaction tx = db.transaction();
				assertThat(tx.numNodesToRead(), is(3L));

				tx.readNode();
				assertThat(tx.state().index, is(0));
				assertThat(tx.conf(), is(conf(-1, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(1024.5)));
				tx.commit();

				assertThat(tx.numNodesToRead(), is(2L));

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(1024.5)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(10.4)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(7.3)));
				assertThat(db.getNumNodes(), is(2L));
			}

			// open the DB and requeue one node
			try (FringeDB db = FringeDB.open(confSpace, file)) {

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(1024.5)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(10.4)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(7.3)));
				assertThat(db.getNumNodes(), is(2L));

				FringeDB.Transaction tx = db.transaction();
				assertThat(tx.numNodesToRead(), is(2L));

				tx.readNode();
				assertThat(tx.state().index, is(1));
				assertThat(tx.conf(), is(conf(-1, -1, -1, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(10.4)));
				tx.writeReplacementNode(tx.state(), tx.conf(), tx.zSumUpper());
				tx.commit();

				assertThat(tx.numNodesToRead(), is(1L));

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(1024.5)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(10.4)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(7.3)));
				assertThat(db.getNumNodes(), is(2L));
			}

			// open the DB and consume one node
			try (FringeDB db = FringeDB.open(confSpace, file)) {

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(1024.5)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(10.4)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(7.3)));
				assertThat(db.getNumNodes(), is(2L));

				FringeDB.Transaction tx = db.transaction();
				assertThat(tx.numNodesToRead(), is(1L));

				tx.readNode();
				assertThat(tx.state().index, is(2));
				assertThat(tx.conf(), is(conf(-1, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(7.3)));
				tx.commit();

				assertThat(tx.numNodesToRead(), is(0L));

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(1024.5)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(10.4)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(7.3)));
				assertThat(db.getNumNodes(), is(1L));
			}

			// open the DB and finish the sweep
			try (FringeDB db = FringeDB.open(confSpace, file)) {

				FringeDB.Transaction tx = db.transaction();
				assertThat(tx.numNodesToRead(), is(0L));
				db.finishStep();

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(Double.NaN)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(10.4)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(Double.NaN)));
				assertThat(db.getNumNodes(), is(1L));
			}

			// check the next sweep
			try (FringeDB db = FringeDB.open(confSpace, file)) {

				FringeDB.Transaction tx = db.transaction();
				assertThat(tx.numNodesToRead(), is(1L));

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(Double.NaN)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(10.4)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(Double.NaN)));
				assertThat(db.getNumNodes(), is(1L));
			}
		}
	}

	@Test
	public void twoLevels() {

		MultiStateConfSpace confSpace = makeConfSpace();

		try (TempFile file = new TempFile("fringe.db")) {
			try (FringeDB db = FringeDB.create(confSpace, file, 1024)) {

				// add the root nodes
				FringeDB.Transaction tx = db.transaction();
				tx.writeRootNode(confSpace.states.get(0), new BigExp(1024.5));
				tx.writeRootNode(confSpace.states.get(1), new BigExp(10.4));
				tx.writeRootNode(confSpace.states.get(2), new BigExp(7.3));
				tx.commit();
				db.finishStep();

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(1024.5)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(10.4)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(7.3)));
				assertThat(db.getNumNodes(), is(3L));

				// do a sweep and replace the roots with children
				tx = db.transaction();

				tx.readNode();
				assertThat(tx.state().index, is(0));
				tx.writeReplacementNode(tx.state(), new int[] { 0, -1 }, new BigExp(35.2));

				tx.readNode();
				assertThat(tx.state().index, is(1));
				tx.writeReplacementNode(tx.state(), new int[] { 1, -1 }, new BigExp(102.3));
				tx.writeReplacementNode(tx.state(), new int[] { -1, 0 }, new BigExp(74.1));

				tx.readNode();
				assertThat(tx.state().index, is(2));

				tx.commit();
				db.finishStep();

				// check db state
				assertThat(db.getZSumMax(confSpace.states.get(0)), is(new BigExp(35.2)));
				assertThat(db.getZSumMax(confSpace.states.get(1)), is(new BigExp(102.3)));
				assertThat(db.getZSumMax(confSpace.states.get(2)), is(new BigExp(Double.NaN)));
				assertThat(db.getNumNodes(), is(3L));

				// sweep the children
				tx = db.transaction();

				tx.readNode();
				assertThat(tx.state().index, is(0));
				assertThat(tx.conf(), is(conf(0, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(35.2)));

				tx.readNode();
				assertThat(tx.state().index, is(1));
				assertThat(tx.conf(), is(conf(1, -1, -1, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(102.3)));

				tx.readNode();
				assertThat(tx.state().index, is(1));
				assertThat(tx.conf(), is(conf(-1, 0, -1, -1)));
				assertThat(tx.zSumUpper(), is(new BigExp(74.1)));

				// TODO: NEXTTIME: make sure the confs agree
			}
		}
	}

	private static MultiStateConfSpace makeConfSpace() {

		Molecule pdb = PDBIO.readResource("/1CC8.ss.pdb");

		Strand design = new Strand.Builder(pdb)
			.setResidues("A68", "A73")
			.build();
		for (String resNum : Arrays.asList("A68", "A69")) {
			design.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}

		Strand target = new Strand.Builder(pdb)
			.setResidues("A2", "A67")
			.build();
		for (String resNum : Arrays.asList("A5", "A6")) {
			target.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}

		// make a multi-state conf space
		return new MultiStateConfSpace
			.Builder("design", new SimpleConfSpace.Builder().addStrands(design).build())
			.addMutableState("complex", new SimpleConfSpace.Builder().addStrands(design, target).build())
			.addUnmutableState("target", new SimpleConfSpace.Builder().addStrands(target).build())
			.build();
	}

	private static Matcher<int[]> conf(int ... expectedRCs) {
		return new BaseMatcher<int[]>() {

			@Override
			public void describeTo(Description desc) {
				desc.appendText("conf")
					.appendValue(Arrays.toString(expectedRCs))
					.appendText("");
			}

			@Override
			public boolean matches(Object item) {
				return item instanceof int[]
					&& Arrays.equals((int[])item, expectedRCs);
			}
		};
	}
}
