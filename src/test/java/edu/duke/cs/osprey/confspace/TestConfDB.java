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

package edu.duke.cs.osprey.confspace;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import com.google.common.collect.Lists;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.TimeTools;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.util.Iterator;
import java.util.function.Consumer;

public class TestConfDB {

	private static SimpleConfSpace confSpace;

	private static File file = new File("conf.db");

	private static SeqSpace.Position lys5;
	private static SeqSpace.Position tyr7;
	private static SeqSpace.Position phe9;

	@BeforeClass
	public static void beforeClass() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A5").setLibraryRotamers(Strand.WildType, "ALA").addWildTypeRotamers();
		strand.flexibility.get("A7").setLibraryRotamers(Strand.WildType, "ALA").addWildTypeRotamers();
		strand.flexibility.get("A9").setLibraryRotamers(Strand.WildType, "ALA").addWildTypeRotamers();

		confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		lys5 = confSpace.seqSpace.getPositionOrThrow("A5");
		tyr7 = confSpace.seqSpace.getPositionOrThrow("A7");
		phe9 = confSpace.seqSpace.getPositionOrThrow("A9");

		assertThat(lys5.wildType.name, is("LYS"));
		assertThat(tyr7.wildType.name, is("TYR"));
		assertThat(phe9.wildType.name, is("PHE"));
	}

	private ConfDB openDB() {
		return new ConfDB(confSpace, file);
	}

	private void cleanDB() {
		if (file.exists()) {
			file.delete();
		}
		assertThat(file.exists(), is(false));
	}

	private void withDB(Consumer<ConfDB> block) {
		cleanDB();
		ConfDB db = openDB();
		try {
			block.accept(db);
		} finally {
			try {
				db.close();
			} catch (Throwable t2) {}
			cleanDB();
		}
	}

	private void withDBTwice(Consumer<ConfDB> block1, Consumer<ConfDB> block2) {
		cleanDB();
		ConfDB db = openDB();
		try {
			block1.accept(db);
			db.close();
		} catch (Throwable t) {
			try {
				db.close();
			} catch (Throwable t2) {}
			cleanDB();
			return;
		}
		db = openDB();
		try {
			block2.accept(db);
		} finally {
			try {
				db.close();
			} catch (Throwable t2) {}
			cleanDB();
		}
	}

	private void assertConf(ConfDB.Conf conf, int[] assignments, double lower, long lowerNs, double upper, long upperNs) {
		assertThat(conf.assignments, is(assignments));
		assertThat(conf.lower.energy, is(lower));
		assertThat(conf.lower.timestampNs, is(lowerNs));
		assertThat(conf.upper.energy, is(upper));
		assertThat(conf.upper.timestampNs, is(upperNs));
	}

	private void assertConfLower(ConfDB.Conf conf, int[] assignments, double lower, long timestampNs) {
		assertThat(conf.assignments, is(assignments));
		assertThat(conf.lower.energy, is(lower));
		assertThat(conf.lower.timestampNs, is(timestampNs));
		assertThat(conf.upper, is(nullValue()));
	}

	private void assertConfUpper(ConfDB.Conf conf, int[] assignments, double upper, long timestampNs) {
		assertThat(conf.assignments, is(assignments));
		assertThat(conf.lower, is(nullValue()));
		assertThat(conf.upper.energy, is(upper));
		assertThat(conf.upper.timestampNs, is(timestampNs));
	}

	@Test
	public void create() {

		cleanDB();

		new ConfDB(confSpace, file).close();

		assertThat(file.exists(), is(true));
		file.delete();
	}

	@Test
	public void createSequenceDB() {
		withDB((db) -> {
			Sequence sequence = confSpace.makeWildTypeSequence();
			ConfDB.SequenceDB sdb = db.getSequence(sequence);
			assertThat(sdb.sequence, sameInstance(sequence));
		});
	}

	@Test
	public void writeReadConfLowerBound() {
		withDB((db) -> {
			Sequence sequence = confSpace.makeWildTypeSequence();
			ConfDB.SequenceDB sdb = db.getSequence(sequence);

			int[] assignments = { 0, 0, 0 };
			double energy = 4.2;
			long timestampNs = TimeTools.getTimestampNs();

			assertThat(assignments.length, is(confSpace.positions.size()));

			sdb.setLowerBound(assignments, energy, timestampNs);

			assertConfLower(sdb.get(assignments), assignments, energy, timestampNs);
		});
	}

	@Test
	public void writeReadConfUpperBound() {
		withDB((db) -> {
			Sequence sequence = confSpace.makeWildTypeSequence();
			ConfDB.SequenceDB sdb = db.getSequence(sequence);

			int[] assignments = { 7, 4, 5 };
			double energy = 7.2;
			long timestampNs = TimeTools.getTimestampNs();

			assertThat(assignments.length, is(confSpace.positions.size()));

			sdb.setUpperBound(assignments, energy, timestampNs);

			assertConfUpper(sdb.get(assignments), assignments, energy, timestampNs);
		});
	}

	@Test
	public void writeReadConfBounds() {
		withDB((db) -> {
			Sequence sequence = confSpace.makeWildTypeSequence();
			ConfDB.SequenceDB sdb = db.getSequence(sequence);

			int[] assignments = { 5, 5, 5 };
			double lowerEnergy = 7.2;
			double upperEnergy = 9.9;
			long timestampNs = TimeTools.getTimestampNs();

			assertThat(assignments.length, is(confSpace.positions.size()));

			sdb.setBounds(assignments, lowerEnergy, upperEnergy, timestampNs);

			assertConf(sdb.get(assignments), assignments, lowerEnergy, timestampNs, upperEnergy, timestampNs);
		});
	}

	@Test
	public void writeReadConfLowerThenUpper() {
		withDB((db) -> {
			Sequence sequence = confSpace.makeWildTypeSequence();
			ConfDB.SequenceDB sdb = db.getSequence(sequence);

			int[] assignments = { 5, 5, 5 };
			double lowerEnergy = 7.2;
			long lowerTimestampNs = TimeTools.getTimestampNs();

			assertThat(assignments.length, is(confSpace.positions.size()));

			sdb.setLowerBound(assignments, lowerEnergy, lowerTimestampNs);

			assertConfLower(sdb.get(assignments), assignments, lowerEnergy, lowerTimestampNs);

			double upperEnergy = 9.9;
			long upperTimestampNs = TimeTools.getTimestampNs();

			sdb.setUpperBound(assignments, upperEnergy, upperTimestampNs);

			assertConf(sdb.get(assignments), assignments, lowerEnergy, lowerTimestampNs, upperEnergy, upperTimestampNs);
		});
	}

	@Test
	public void writeReadAFewConfs() {
		withDB((db) -> {
			Sequence sequence = confSpace.makeWildTypeSequence();
			ConfDB.SequenceDB sdb = db.getSequence(sequence);

			sdb.setUpperBound(new int[] { 1, 2, 3 }, 7.9, 42L);
			sdb.setUpperBound(new int[] { 7, 9, 8 }, 3.2, 54L);
			sdb.setUpperBound(new int[] { 4, 0, 5 }, 2.3, 69L);

			Iterator<ConfDB.Conf> confs = sdb.iterator();

			// confs should come out in lexicographic order of the assignments
			assertConfUpper(confs.next(), new int[] { 1, 2, 3 }, 7.9, 42L);
			assertConfUpper(confs.next(), new int[] { 4, 0, 5 }, 2.3, 69L);
			assertConfUpper(confs.next(), new int[] { 7, 9, 8 }, 3.2, 54L);
			assertThat(confs.hasNext(), is(false));
		});
	}

	@Test
	public void writeCloseReadAFewConfs() {
		Sequence sequence = confSpace.makeWildTypeSequence();
		withDBTwice((db) -> {

			ConfDB.SequenceDB sdb = db.getSequence(sequence);
			sdb.setUpperBound(new int[] { 1, 2, 3 }, 7.9, 42L);
			sdb.setUpperBound(new int[] { 7, 9, 8 }, 3.2, 54L);
			sdb.setUpperBound(new int[] { 4, 0, 5 }, 2.3, 69L);

		}, (db) -> {

			Iterator<ConfDB.Conf> confs = db.getSequence(sequence).iterator();

			// confs should come out in lexicographic order of the assignments
			assertConfUpper(confs.next(), new int[] { 1, 2, 3 }, 7.9, 42L);
			assertConfUpper(confs.next(), new int[] { 4, 0, 5 }, 2.3, 69L);
			assertConfUpper(confs.next(), new int[] { 7, 9, 8 }, 3.2, 54L);
			assertThat(confs.hasNext(), is(false));
		});
	}

	@Test
	public void writeReadCloseReadAFewSequences() {

		Sequence sequence1 = confSpace.makeUnassignedSequence()
			.set(lys5, "LYS")
			.set(tyr7, "TYR")
			.set(phe9, "PHE");

		Sequence sequence2 = confSpace.makeUnassignedSequence()
			.set(lys5, "ALA")
			.set(tyr7, "TYR")
			.set(phe9, "PHE");

		Sequence sequence3 = confSpace.makeUnassignedSequence()
			.set(lys5, "LYS")
			.set(tyr7, "ALA")
			.set(phe9, "PHE");

		Sequence sequence4 = confSpace.makeUnassignedSequence()
			.set(lys5, "LYS")
			.set(tyr7, "TYR")
			.set(phe9, "ALA");

		withDBTwice((db) -> {

			ConfDB.SequenceDB sdb1 = db.getSequence(sequence1);
			ConfDB.SequenceDB sdb2 = db.getSequence(sequence2);
			ConfDB.SequenceDB sdb3 = db.getSequence(sequence3);
			ConfDB.SequenceDB sdb4 = db.getSequence(sequence4);

			assertThat(sdb1.sequence, sameInstance(sequence1));
			assertThat(sdb2.sequence, sameInstance(sequence2));
			assertThat(sdb3.sequence, sameInstance(sequence3));
			assertThat(sdb4.sequence, sameInstance(sequence4));

			// sequences get hased in the db, so they can come out in any order
			assertThat(Lists.newArrayList(db.getSequences()), containsInAnyOrder(
				sequence1, sequence2, sequence3, sequence4
			));

		}, (db) -> {

			assertThat(Lists.newArrayList(db.getSequences()), containsInAnyOrder(
				sequence1, sequence2, sequence3, sequence4
			));
		});
	}

	@Test
	public void writeCloseReadSequenceInfo() {

		Sequence sequence = confSpace.makeWildTypeSequence();

		withDBTwice((db) -> {

			ConfDB.SequenceDB sdb = db.getSequence(sequence);
			assertThat(sdb.getLowerEnergyOfUnsampledConfs(), is(Double.NaN));

			sdb.setLowerEnergyOfUnsampledConfs(4.2);

			assertThat(sdb.getLowerEnergyOfUnsampledConfs(), is(4.2));

		}, (db) -> {

			ConfDB.SequenceDB sdb = db.getSequence(sequence);
			assertThat(sdb.getLowerEnergyOfUnsampledConfs(), is(4.2));
		});
	}

	@Test
	public void astarSearch() {
		withDBTwice((db) -> {

			new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
				.setParallelism(Parallelism.makeCpu(4))
				.use((ecalc) -> {

					ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

					EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
						.build()
						.calcEnergyMatrix();

					ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
						.setTraditional()
						.setShowProgress(true)
						.build();


					long numConfs = astar.getNumConformations().longValueExact();
					assertThat(numConfs, is(1740L));

					// write all the conformations in the A* tree
					ConfSearch.ScoredConf conf;
					while ((conf = astar.nextConf()) != null) {
						ConfDB.SequenceDB sdb = db.getSequence(confSpace.makeSequenceFromConf(conf));
						sdb.setLowerBound(
							conf.getAssignments(),
							conf.getScore(),
							TimeTools.getTimestampNs()
						);
						sdb.updateLowerEnergyOfUnsampledConfs(conf.getScore());

					}
				});

		}, (db) -> {

			assertThat(db.getNumSequences(), is(8L)); // 2^3 sequences

			// count the confs
			long numConfs = 0;
			for (Sequence sequence : db.getSequences()) {
				for (ConfDB.Conf conf : db.getSequence(sequence)) {
					numConfs++;
				}
			}
			assertThat(numConfs, is(1740L));

			// check the lower bounds
			for (Sequence sequence : db.getSequences()) {
				ConfDB.SequenceDB sdb = db.getSequence(sequence);
				double lowerEnergy = Double.NEGATIVE_INFINITY;
				for (ConfDB.Conf conf : sdb) {
					lowerEnergy = Math.max(lowerEnergy, conf.lower.energy);
				}
				assertThat(lowerEnergy, is(sdb.getLowerEnergyOfUnsampledConfs()));
			}
		});
	}

	@Test
	public void energyIndices() {

		int[][] assignments = {
			{ 0, 0, 0 },
			{ 1, 2, 3 },
			{ 3, 2, 1 },
		};

		String tableId = "foo";
		withDBTwice((db) -> {

			ConfDB.ConfTable table = db.new ConfTable(tableId);

			table.setBounds(assignments[0], 7.0, 27.0, 5L);
			table.setBounds(assignments[1], 6.0, 25.0, 6L);
			table.setBounds(assignments[2], 5.0, 26.0, 7L);

			assertThat(table.size(), is(3L));

			assertThat(table.energiedConfs(ConfDB.SortOrder.Assignment), contains(
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 27.0),
				new ConfSearch.EnergiedConf(assignments[1], 6.0, 25.0),
				new ConfSearch.EnergiedConf(assignments[2], 5.0, 26.0)
			));

			assertThat(table.energiedConfs(ConfDB.SortOrder.Score), contains(
				new ConfSearch.EnergiedConf(assignments[2], 5.0, 26.0),
				new ConfSearch.EnergiedConf(assignments[1], 6.0, 25.0),
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 27.0)
			));

			assertThat(table.energiedConfs(ConfDB.SortOrder.Energy), contains(
				new ConfSearch.EnergiedConf(assignments[1], 6.0, 25.0),
				new ConfSearch.EnergiedConf(assignments[2], 5.0, 26.0),
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 27.0)
			));

			assertThat(table.lowerBounds(), contains(5.0, 6.0, 7.0));
			assertThat(table.upperBounds(), contains(25.0, 26.0, 27.0));

			Iterator<ConfDB.Conf> iter;

			iter = table.getConfsByLowerBound(5.0).iterator();
			assertConf(iter.next(), assignments[2], 5.0, 7L, 26.0, 7L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByLowerBound(6.0).iterator();
			assertConf(iter.next(), assignments[1], 6.0, 6L, 25.0, 6L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByLowerBound(7.0).iterator();
			assertConf(iter.next(), assignments[0], 7.0, 5L, 27.0, 5L);
			assertThat(iter.hasNext(), is(false));

			assertThat(table.getConfsByLowerBound(4.0), is(nullValue()));
			assertThat(table.getConfsByLowerBound(8.0), is(nullValue()));

			iter = table.getConfsByUpperBound(25.0).iterator();
			assertConf(iter.next(), assignments[1], 6.0, 6L, 25.0, 6L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByUpperBound(26.0).iterator();
			assertConf(iter.next(), assignments[0], 7.0, 5L, 27.0, 5L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByUpperBound(27.0).iterator();
			assertConf(iter.next(), assignments[2], 5.0, 7L, 26.0, 7L);
			assertThat(iter.hasNext(), is(false));

			assertThat(table.getConfsByUpperBound(24.0), is(nullValue()));
			assertThat(table.getConfsByUpperBound(28.0), is(nullValue()));

		}, (db) -> {

			ConfDB.ConfTable table = db.new ConfTable(tableId);

			assertThat(table.size(), is(3L));

			assertThat(table.energiedConfs(ConfDB.SortOrder.Assignment), contains(
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 27.0),
				new ConfSearch.EnergiedConf(assignments[1], 6.0, 25.0),
				new ConfSearch.EnergiedConf(assignments[2], 5.0, 26.0)
			));

			assertThat(table.energiedConfs(ConfDB.SortOrder.Score), contains(
				new ConfSearch.EnergiedConf(assignments[2], 5.0, 26.0),
				new ConfSearch.EnergiedConf(assignments[1], 6.0, 25.0),
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 27.0)
			));

			assertThat(table.energiedConfs(ConfDB.SortOrder.Energy), contains(
				new ConfSearch.EnergiedConf(assignments[1], 6.0, 25.0),
				new ConfSearch.EnergiedConf(assignments[2], 5.0, 26.0),
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 27.0)
			));

			assertThat(table.lowerBounds(), contains(5.0, 6.0, 7.0));
			assertThat(table.upperBounds(), contains(25.0, 26.0, 27.0));

			Iterator<ConfDB.Conf> iter;

			iter = table.getConfsByLowerBound(5.0).iterator();
			assertConf(iter.next(), assignments[2], 5.0, 7L, 26.0, 7L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByLowerBound(6.0).iterator();
			assertConf(iter.next(), assignments[1], 6.0, 6L, 25.0, 6L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByLowerBound(7.0).iterator();
			assertConf(iter.next(), assignments[0], 7.0, 5L, 27.0, 5L);
			assertThat(iter.hasNext(), is(false));

			assertThat(table.getConfsByLowerBound(4.0), is(nullValue()));
			assertThat(table.getConfsByLowerBound(8.0), is(nullValue()));

			iter = table.getConfsByUpperBound(25.0).iterator();
			assertConf(iter.next(), assignments[1], 6.0, 6L, 25.0, 6L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByUpperBound(26.0).iterator();
			assertConf(iter.next(), assignments[0], 7.0, 5L, 27.0, 5L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByUpperBound(27.0).iterator();
			assertConf(iter.next(), assignments[2], 5.0, 7L, 26.0, 7L);
			assertThat(iter.hasNext(), is(false));

			assertThat(table.getConfsByUpperBound(24.0), is(nullValue()));
			assertThat(table.getConfsByUpperBound(28.0), is(nullValue()));
		});
	}

	@Test
	public void energyIndicesChange() {

		int[][] assignments = {
			{ 0, 0, 0 },
			{ 1, 2, 3 },
			{ 3, 2, 1 },
		};

		withDB((db) -> {

			ConfDB.ConfTable table = db.new ConfTable("foo");

			table.setBounds(assignments[0], 7.0, 27.0, 5L);
			table.setBounds(assignments[1], 6.0, 25.0, 6L);
			table.setBounds(assignments[2], 5.0, 26.0, 7L);

			table.setBounds(assignments[1], 10.0, 50.0, 10L);

			assertThat(table.size(), is(3L));

			assertThat(table.energiedConfs(ConfDB.SortOrder.Assignment), contains(
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 27.0),
				new ConfSearch.EnergiedConf(assignments[1], 10.0, 50.0),
				new ConfSearch.EnergiedConf(assignments[2], 5.0, 26.0)
			));

			assertThat(table.energiedConfs(ConfDB.SortOrder.Score), contains(
				new ConfSearch.EnergiedConf(assignments[2], 5.0, 26.0),
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 27.0),
				new ConfSearch.EnergiedConf(assignments[1], 10.0, 50.0)
			));

			assertThat(table.energiedConfs(ConfDB.SortOrder.Energy), contains(
				new ConfSearch.EnergiedConf(assignments[2], 5.0, 26.0),
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 27.0),
				new ConfSearch.EnergiedConf(assignments[1], 10.0, 50.0)
			));

			assertThat(table.lowerBounds(), contains(5.0, 7.0, 10.0));
			assertThat(table.upperBounds(), contains(26.0, 27.0, 50.0));

			Iterator<ConfDB.Conf> iter;

			iter = table.getConfsByLowerBound(5.0).iterator();
			assertConf(iter.next(), assignments[2], 5.0, 7L, 26.0, 7L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByLowerBound(7.0).iterator();
			assertConf(iter.next(), assignments[0], 7.0, 5L, 27.0, 5L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByLowerBound(10.0).iterator();
			assertConf(iter.next(), assignments[1], 10.0, 10L, 50.0, 10L);
			assertThat(iter.hasNext(), is(false));

			assertThat(table.getConfsByLowerBound(6.0), is(nullValue()));

			iter = table.getConfsByUpperBound(26.0).iterator();
			assertConf(iter.next(), assignments[2], 5.0, 7L, 26.0, 7L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByUpperBound(27.0).iterator();
			assertConf(iter.next(), assignments[0], 7.0, 5L, 27.0, 5L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByUpperBound(50.0).iterator();
			assertConf(iter.next(), assignments[1], 10.0, 10L, 50.0, 10L);
			assertThat(iter.hasNext(), is(false));

			assertThat(table.getConfsByUpperBound(25.0), is(nullValue()));
		});
	}

	@Test
	public void energyIndexChangeLowerBound() {

		int[][] assignments = {
			{ 0, 0, 0 },
			{ 1, 2, 3 },
			{ 3, 2, 1 },
		};

		withDB((db) -> {

			ConfDB.ConfTable table = db.new ConfTable("foo");

			table.setLowerBound(assignments[0], 7.0, 5L);
			table.setLowerBound(assignments[1], 6.0, 6L);
			table.setLowerBound(assignments[2], 5.0, 7L);

			table.setLowerBound(assignments[1], 10.0, 10L);

			assertThat(table.scoredConfs(ConfDB.SortOrder.Assignment), contains(
				new ConfSearch.ScoredConf(assignments[0], 7.0),
				new ConfSearch.ScoredConf(assignments[1], 10.0),
				new ConfSearch.ScoredConf(assignments[2], 5.0)
			));

			assertThat(table.scoredConfs(ConfDB.SortOrder.Score), contains(
				new ConfSearch.ScoredConf(assignments[2], 5.0),
				new ConfSearch.ScoredConf(assignments[0], 7.0),
				new ConfSearch.ScoredConf(assignments[1], 10.0)
			));

			assertThat(table.lowerBounds(), contains(5.0, 7.0, 10.0));

			Iterator<ConfDB.Conf> iter;

			iter = table.getConfsByLowerBound(5.0).iterator();
			assertConfLower(iter.next(), assignments[2], 5.0, 7L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByLowerBound(7.0).iterator();
			assertConfLower(iter.next(), assignments[0], 7.0, 5L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByLowerBound(10.0).iterator();
			assertConfLower(iter.next(), assignments[1], 10.0, 10L);
			assertThat(iter.hasNext(), is(false));

			assertThat(table.getConfsByLowerBound(6.0), is(nullValue()));
		});
	}

	@Test
	public void energyIndexChangeUpperBound() {

		int[][] assignments = {
			{ 0, 0, 0 },
			{ 1, 2, 3 },
			{ 3, 2, 1 },
		};

		withDB((db) -> {

			ConfDB.ConfTable table = db.new ConfTable("foo");

			table.setUpperBound(assignments[0], 27.0, 5L);
			table.setUpperBound(assignments[1], 25.0, 6L);
			table.setUpperBound(assignments[2], 26.0, 7L);

			table.setUpperBound(assignments[1], 50.0, 10L);

			assertThat(table.energiedConfs(ConfDB.SortOrder.Assignment), contains(
				new ConfSearch.EnergiedConf(assignments[0], Double.NaN, 27.0),
				new ConfSearch.EnergiedConf(assignments[1], Double.NaN, 50.0),
				new ConfSearch.EnergiedConf(assignments[2], Double.NaN, 26.0)
			));

			assertThat(table.energiedConfs(ConfDB.SortOrder.Energy), contains(
				new ConfSearch.EnergiedConf(assignments[2], Double.NaN, 26.0),
				new ConfSearch.EnergiedConf(assignments[0], Double.NaN, 27.0),
				new ConfSearch.EnergiedConf(assignments[1], Double.NaN, 50.0)
			));

			assertThat(table.upperBounds(), contains(26.0, 27.0, 50.0));

			Iterator<ConfDB.Conf> iter;

			iter = table.getConfsByUpperBound(26.0).iterator();
			assertConfUpper(iter.next(), assignments[2], 26.0, 7L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByUpperBound(27.0).iterator();
			assertConfUpper(iter.next(), assignments[0], 27.0, 5L);
			assertThat(iter.hasNext(), is(false));

			iter = table.getConfsByUpperBound(50.0).iterator();
			assertConfUpper(iter.next(), assignments[1], 50.0, 10L);
			assertThat(iter.hasNext(), is(false));

			assertThat(table.getConfsByUpperBound(25.0), is(nullValue()));
		});
	}

	@Test
	public void energyIndicesSameEnergies() {

		int[][] assignments = {
			{ 0, 0, 0 },
			{ 1, 2, 3 },
			{ 3, 2, 1 },
		};

		String tableId = "foo";
		withDBTwice((db) -> {

			ConfDB.ConfTable table = db.new ConfTable(tableId);

			table.setBounds(assignments[0], 7.0, 20.0, 5L);
			table.setBounds(assignments[1], 7.0, 20.0, 6L);
			table.setBounds(assignments[2], 7.0, 20.0, 7L);

			assertThat(table.size(), is(3L));

			assertThat(table.scoredConfs(ConfDB.SortOrder.Assignment), containsInAnyOrder(
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 20.0),
				new ConfSearch.EnergiedConf(assignments[1], 7.0, 20.0),
				new ConfSearch.EnergiedConf(assignments[2], 7.0, 20.0)
			));

			assertThat(table.scoredConfs(ConfDB.SortOrder.Score), containsInAnyOrder(
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 20.0),
				new ConfSearch.EnergiedConf(assignments[1], 7.0, 20.0),
				new ConfSearch.EnergiedConf(assignments[2], 7.0, 20.0)
			));

			assertThat(table.scoredConfs(ConfDB.SortOrder.Energy), containsInAnyOrder(
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 20.0),
				new ConfSearch.EnergiedConf(assignments[1], 7.0, 20.0),
				new ConfSearch.EnergiedConf(assignments[2], 7.0, 20.0)
			));

			assertThat(table.lowerBounds(), contains(7.0));
			assertThat(table.lowerBounds(), contains(20.0));

			// order should be add order
			Iterator<ConfDB.Conf> iter = table.getConfsByLowerBound(7.0).iterator();
			assertConf(iter.next(), assignments[0], 7.0, 5L, 20.0, 5L);
			assertConf(iter.next(), assignments[1], 7.0, 6L, 20.0, 6L);
			assertConf(iter.next(), assignments[2], 7.0, 7L, 20.0, 7L);
			assertThat(iter.hasNext(), is(false));

			assertThat(table.getConfsByLowerBound(5.0).iterator().hasNext(), is(false));
			assertThat(table.getConfsByLowerBound(6.0).iterator().hasNext(), is(false));

		}, (db) -> {

			ConfDB.ConfTable table = db.new ConfTable(tableId);

			assertThat(table.size(), is(3L));

			assertThat(table.scoredConfs(ConfDB.SortOrder.Assignment), containsInAnyOrder(
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 20.0),
				new ConfSearch.EnergiedConf(assignments[1], 7.0, 20.0),
				new ConfSearch.EnergiedConf(assignments[2], 7.0, 20.0)
			));

			assertThat(table.scoredConfs(ConfDB.SortOrder.Score), containsInAnyOrder(
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 20.0),
				new ConfSearch.EnergiedConf(assignments[1], 7.0, 20.0),
				new ConfSearch.EnergiedConf(assignments[2], 7.0, 20.0)
			));

			assertThat(table.scoredConfs(ConfDB.SortOrder.Energy), containsInAnyOrder(
				new ConfSearch.EnergiedConf(assignments[0], 7.0, 20.0),
				new ConfSearch.EnergiedConf(assignments[1], 7.0, 20.0),
				new ConfSearch.EnergiedConf(assignments[2], 7.0, 20.0)
			));

			// order should be add order
			Iterator<ConfDB.Conf> iter = table.getConfsByLowerBound(7.0).iterator();
			assertConf(iter.next(), assignments[0], 7.0, 5L, 20.0, 5L);
			assertConf(iter.next(), assignments[1], 7.0, 6L, 20.0, 6L);
			assertConf(iter.next(), assignments[2], 7.0, 7L, 20.0, 7L);
			assertThat(iter.hasNext(), is(false));

			assertThat(table.getConfsByLowerBound(5.0).iterator().hasNext(), is(false));
			assertThat(table.getConfsByLowerBound(6.0).iterator().hasNext(), is(false));
		});
	}
}
