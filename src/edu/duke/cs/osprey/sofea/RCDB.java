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

import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.tools.*;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;

import java.io.File;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Consumer;


public class RCDB implements AutoCloseable {

	public final MultiStateConfSpace confSpace;
	public final MathContext mathContext;
	public final File file;

	private final DB db;
	private final Map<String,Table> tables = new HashMap<>();

	public RCDB(MultiStateConfSpace confSpace, MathContext mathContext) {
		this(confSpace, mathContext, null);
	}

	public RCDB(MultiStateConfSpace confSpace, MathContext mathContext, File file) {

		this.confSpace = confSpace;
		this.mathContext = mathContext;
		this.file = file;

		// open the DB
		if (file != null) {
			db = DBMaker.fileDB(file)
				.fileMmapEnableIfSupported() // use memory-mapped files if possible (can be much faster)
				.closeOnJvmShutdown()
				.make();
		} else {
			db = DBMaker.memoryDB()
				.make();
		}
	}

	@Override
	public void close() {
		db.close();
	}

	public BigMath bigMath() {
		return new BigMath(mathContext);
	}

	public class Table {

		public final MultiStateConfSpace.State state;
		public final Sequence seq;
		public final String id;

		private final int[][] rcIds;
		private final HTreeMap<Integer,BigDecimalBounds> map;

		public Table(MultiStateConfSpace.State state, Sequence seq, String id) {

			this.state = state;
			this.seq = seq;
			this.id = id;

			// assign ids to all the rcs
			rcIds = new int[state.confSpace.positions.size()][];
			int nextId = 0;
			for (SimpleConfSpace.Position pos : state.confSpace.positions) {
				rcIds[pos.index] = new int[pos.resConfs.size()];
				for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {
					rcIds[pos.index][rc.index] = nextId++;
				}
			}

			// open the table
			map = db.hashMap(id)
				.keySerializer(Serializer.INTEGER)
				.valueSerializer(new MapDBTools.BigDecimalBoundsSerializer())
				.createOrOpen();
		}

		private int getId(int pos, int rc) {
			return rcIds[pos][rc];
		}

		public BigDecimalBounds getZSumBounds(int pos, int rc) {
			// get the sum, or [0,0]
			BigDecimalBounds sum = map.get(getId(pos, rc));
			if (sum == null) {
				sum = new BigDecimalBounds(BigDecimal.ZERO, BigDecimal.ZERO);
			}
			return sum;
		}

		private boolean isZero(BigDecimalBounds bounds) {
			return MathTools.isZero(bounds.lower)
				&& MathTools.isZero(bounds.upper);
		}

		public void update(int pos, int rc, Consumer<BigDecimalBounds> updater) {

			BigDecimalBounds sum = getZSumBounds(pos, rc);

			updater.accept(sum);

			// save the sum, if it's not [0,0]
			if (isZero(sum)) {
				map.remove(getId(pos, rc));
			} else {
				map.put(getId(pos, rc), sum);
			}
		}

		public void addZSumUpper(int pos, int rc, BigDecimal zSumUpper) {

			if (!MathTools.isFinite(zSumUpper)) {
				throw new IllegalArgumentException("Z must be finite: " + zSumUpper);
			}

			update(pos, rc, sum -> {
				sum.upper = bigMath()
					.set(sum.upper)
					.add(zSumUpper)
					.get();
			});
		}

		public void subZSumUpper(int pos, int rc, BigDecimal zSumUpper) {

			if (!MathTools.isFinite(zSumUpper)) {
				throw new IllegalArgumentException("Z must be finite: " + zSumUpper);
			}

			update(pos, rc, sum -> {
				sum.upper = bigMath()
					.set(sum.upper)
					.sub(zSumUpper)
					.get();
			});
		}

		public void addZPath(int[] conf, int pos, BigDecimal zPath, BigDecimal zSumUpper) {

			if (!MathTools.isFinite(zPath) || !MathTools.isFinite(zSumUpper)) {
				throw new IllegalArgumentException("Z must be finite: " + zPath + ", " + zSumUpper);
			}

			for (SimpleConfSpace.Position opos : state.confSpace.positions) {
				int rc = conf[opos.index];

				update(opos.index, rc, sum -> {
					sum.lower = bigMath()
						.set(sum.lower)
						.add(zPath)
						.get();
					sum.upper = bigMath()
						.set(sum.upper)
						.add(zPath)
						.get();
				});
			}

			// subtract off the zSumUpper only on the last position to be assigned (ie, the leaf node)
			update(pos, conf[pos], sum ->
				sum.upper = bigMath()
					.set(sum.upper)
					.sub(zSumUpper)
					.get()
			);
		}
	}

	public Table table(MultiStateConfSpace.State state, Sequence seq) {
		String id = String.format("%d[%s]", state.index, seq.toString(Sequence.Renderer.ResType));
		return tables.computeIfAbsent(id, key -> new Table(state, seq, id));
	}

	public void addZSumUpper(MultiStateConfSpace.State state, Sequence seq, int[] conf, int pos, BigDecimal zSumUpper) {
		table(state, seq).addZSumUpper(pos, conf[pos], zSumUpper);
	}

	public void subZSumUpper(MultiStateConfSpace.State state, Sequence seq, int[] conf, int pos, BigDecimal zSumUpper) {
		table(state, seq).subZSumUpper(pos, conf[pos], zSumUpper);
	}

	public void addZPath(MultiStateConfSpace.State state, Sequence seq, int[] conf, int pos, BigDecimal zPath, BigDecimal zSumUpper) {
		table(state, seq).addZPath(conf, pos, zPath, zSumUpper);
	}

	public BigDecimalBounds getZSumBounds(MultiStateConfSpace.State state, Sequence seq, int pos, int rc) {

		// start with the bounds from the full sequence at this pos,rc
		BigDecimalBounds zSumBounds = table(state, seq).getZSumBounds(pos, rc);
		BigMath upper = bigMath().set(zSumBounds.upper);

		// (the lower bounds are zero)
		forEachPartialSequence(seq, partialSeq -> {
			BigDecimalBounds partialZSumBounds = table(state, partialSeq).getZSumBounds(pos, rc);
			upper.add(partialZSumBounds.upper);
		});

		return new BigDecimalBounds(
			zSumBounds.lower,
			upper.get()
		);
	}

	private void forEachPartialSequence(Sequence seq, Consumer<Sequence> block) {

		// iterate backwards and un-assign residues
		// NOTE: assumes residues are assigned in seq space position order
		Sequence parentSeq = seq.copy();
		for (int i=seq.countAssignments() - 1; i>=0; i--) {
			parentSeq.rtIndices[i] = Sequence.Unassigned;
			block.accept(parentSeq);
		}
	}

	// TODO: implement transactions
	public void commit() {
		db.commit();
	}
}
