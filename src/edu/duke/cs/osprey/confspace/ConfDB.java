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

import edu.duke.cs.osprey.tools.AutoCleanable;
import edu.duke.cs.osprey.tools.IntEncoding;
import edu.duke.cs.osprey.tools.Streams;

import edu.duke.cs.osprey.tools.UnpossibleError;
import org.jetbrains.annotations.NotNull;
import org.mapdb.*;
import org.mapdb.serializer.GroupSerializerObjectArray;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;


public class ConfDB implements AutoCleanable {

	public static ConfDB makeIfNeeded(SimpleConfSpace confSpace, File file) {

		// no file? confdb not needed
		if (file == null) {
			return null;
		}

		return new ConfDB(confSpace, file);
	}

	public static class DBs implements AutoCleanable {

		public class Adder {

			public void add(SimpleConfSpace confSpace, File file) {
				if (file != null) {
					dbs.put(confSpace, new ConfDB(confSpace, file));
				}
			}
		}

		private final Map<SimpleConfSpace,ConfDB> dbs = new HashMap<>();
		private final Adder adder = new Adder();

		public DBs add(SimpleConfSpace confSpace, File file) {
			adder.add(confSpace, file);
			return this;
		}

		public <T> DBs addAll(Iterable<T> things, BiConsumer<T,Adder> block) {
			for (T thing : things) {
				block.accept(thing, adder);
			}
			return this;
		}

		public ConfDB get(SimpleConfSpace confSpace) {
			return dbs.get(confSpace);
		}

		@Override
		public void clean() {
			for (ConfDB db : dbs.values()) {
				db.clean();
			}
		}
	}

	public static enum SortOrder {
		Assignment,
		Score,
		Energy
	}

	public static class Conf {

		public static class Bound {

			public final double energy;
			public final long timestampNs;

			public Bound(double energy, long timestampNs) {
				this.energy = energy;
				this.timestampNs = timestampNs;
			}
		}

		public final int[] assignments;
		public final Bound lower;
		public final Bound upper;

		public Conf(int[] assignments, Bound lower, Bound upper) {
			this.assignments = assignments;
			this.lower = lower;
			this.upper = upper;
		}

		private Conf(int[] assignments, ConfInfo info) {
			this.assignments = assignments;
			this.lower = info.makeLowerBound();
			this.upper = info.makeUpperBound();
		}

		public ConfSearch.ScoredConf toScoredConf() {
			return new ConfSearch.ScoredConf(
				assignments,
				lower == null ? Double.NaN : lower.energy
			);
		}

		public ConfSearch.EnergiedConf toEnergiedConf() {
			return new ConfSearch.EnergiedConf(
				assignments,
				lower == null ? Double.NaN : lower.energy,
				upper == null ? Double.NaN : upper.energy
			);
		}
	}

	private static class ConfInfo {

		public double lowerEnergy;
		public long lowerTimestampNs;
		public double upperEnergy;
		public long upperTimestampNs;

		public ConfInfo() {
			lowerEnergy = 0.0;
			lowerTimestampNs = 0L;
			upperEnergy = 0.0;
			upperTimestampNs = 0L;
		}

		public ConfInfo(double lowerEnergy, long lowerTimestampNs, double upperEnergy, long upperTimestampNs) {
			this.lowerEnergy = lowerEnergy;
			this.lowerTimestampNs = lowerTimestampNs;
			this.upperEnergy = upperEnergy;
			this.upperTimestampNs = upperTimestampNs;
		}

		public Conf.Bound makeLowerBound() {
			return makeBound(lowerEnergy, lowerTimestampNs);
		}

		public Conf.Bound makeUpperBound() {
			return makeBound(upperEnergy, upperTimestampNs);
		}

		private Conf.Bound makeBound(double energy, long timestampNs) {
			if (timestampNs == 0L) {
				return null;
			} else {
				return new Conf.Bound(energy, timestampNs);
			}
		}
	}

	private static abstract class SimpleSerializer<T> extends GroupSerializerObjectArray<T> {

		public static final int DynamicSize = -1;

		public final int fixedSize;

		protected SimpleSerializer(int fixedSize) {
			this.fixedSize = fixedSize;
		}

		@Override
		public boolean isTrusted() {
			// we always read/write the same number of bytes, so we're "trusted" by MapDB
			return true;
		}

		@Override
		public int fixedSize() {
			return fixedSize;
		}

		@Override
		public T deserializeFromLong(long input, int available) {
			// nope
			throw new UnsupportedOperationException();
		}
	}

	private class AssignmentsSerializer extends SimpleSerializer<int[]> {

		private final int numPos;

		public AssignmentsSerializer() {
			super(assignmentEncoding.numBytes*confSpace.positions.size());
			this.numPos = confSpace.positions.size();

			// check the unassigned value is -1, since we do arithmetic on it
			assert (edu.duke.cs.osprey.confspace.Conf.Unassigned == -1);
		}

		@Override
		public void serialize(@NotNull DataOutput2 out, @NotNull int[] assignments)
		throws IOException {
			for (int i=0; i<numPos; i++) {
				assignmentEncoding.write(out, assignments[i] + 1); // +1 to shift the unassigned value (-1) to non-negative
			}
		}

		@Override
		public int[] deserialize(@NotNull DataInput2 in, int available)
		throws IOException {
			int[] assignments = new int[numPos];
			for (int i=0; i<numPos; i++) {
				assignments[i] = assignmentEncoding.read(in) - 1;
			}
			return assignments;
		}

		@Override
		public int compare(int[] a, int[] b) {

			// short circuit
			if (a == b) {
				return 0;
			}

			// lexicographical comparison
			final int len = Math.min(a.length, b.length);
			for (int i=0; i<len; i++) {
				int val = Integer.compare(a[i], b[i]);
				if (val != 0) {
					return val;
				}
			}
			return Integer.compare(a.length, b.length);
		}

		@Override
		public boolean equals(int[] a, int[] b) {
			return Arrays.equals(a, b);
		}

		@Override
		public int hashCode(@NotNull int[] assignments, int seed) {
			return DataIO.intHash(Arrays.hashCode(assignments) + seed);
		}
	}

	private class MultiAssignmentsSerializer extends SimpleSerializer<List<int[]>> {

		private final AssignmentsSerializer serializer = new AssignmentsSerializer();

		public MultiAssignmentsSerializer() {
			super(SimpleSerializer.DynamicSize);
		}

		@Override
		public void serialize(@NotNull DataOutput2 out, @NotNull List<int[]> multiAssignments)
		throws IOException {
			out.writeInt(multiAssignments.size());
			for (int[] assignments : multiAssignments) {
				serializer.serialize(out, assignments);
			}
		}

		@Override
		public List<int[]> deserialize(@NotNull DataInput2 in, int available)
		throws IOException {
			List<int[]> multiAssignments = new ArrayList<>();
			int num = in.readInt();
			for (int i=0; i<num; i++) {
				multiAssignments.add(serializer.deserialize(in, available));
			}
			return multiAssignments;
		}
	}

	private static class SequenceInfo {

		public double lowerEnergyOfUnsampledConfs;
		// TODO: other sequence-level properties?

		public SequenceInfo() {
			this.lowerEnergyOfUnsampledConfs = Double.NaN;
		}

		public SequenceInfo(double lowerEnergyOfUnsampledConfs) {
			this.lowerEnergyOfUnsampledConfs = lowerEnergyOfUnsampledConfs;
		}
	}

	public class ConfTable implements Iterable<Conf> {

		private final BTreeMap<int[],ConfInfo> btree;
		private final EnergyIndex lowerIndex;
		private final EnergyIndex upperIndex;

		public ConfTable(String id) {

			// MapDB serializer for ConfInfo
			final int ConfInfoBytes = Double.BYTES*2 + Long.BYTES*2;
			SimpleSerializer<ConfInfo> confInfoSerializer = new SimpleSerializer<ConfInfo>(ConfInfoBytes) {

				@Override
				public void serialize(@NotNull DataOutput2 out, @NotNull ConfInfo info)
				throws IOException {
					out.writeDouble(info.lowerEnergy);
					out.writeLong(info.lowerTimestampNs);
					out.writeDouble(info.upperEnergy);
					out.writeLong(info.upperTimestampNs);
				}

				@Override
				public ConfInfo deserialize(@NotNull DataInput2 in, int available)
				throws IOException {
					return new ConfInfo(
						in.readDouble(),
						in.readLong(),
						in.readDouble(),
						in.readLong()
					);
				}
			};

			this.btree = db.treeMap(id)
				.keySerializer(new AssignmentsSerializer())
				.valueSerializer(confInfoSerializer)
				.createOrOpen();

			this.lowerIndex = new EnergyIndex(id + "-lowerEnergy");
			this.upperIndex = new EnergyIndex(id + "-upperEnergy");
		}

		public void setBounds(ConfSearch.EnergiedConf econf, long timestampNs) {
			setBounds(econf.getAssignments(), econf.getScore(), econf.getEnergy(), timestampNs);
		}

		public void setBounds(int[] assignments, double lowerEnergy, double upperEnergy, long timestampNs) {
			ConfInfo info = btree.get(assignments);
			if (info == null) {
				info = new ConfInfo();
			} else {
				// remove old energy index entries if needed
				if (info.lowerTimestampNs != 0L) {
					lowerIndex.remove(info.lowerEnergy, assignments);
				}
				if (info.upperTimestampNs != 0L) {
					upperIndex.remove(info.upperEnergy, assignments);
				}
			}
			info.lowerEnergy = lowerEnergy;
			info.lowerTimestampNs = timestampNs;
			info.upperEnergy = upperEnergy;
			info.upperTimestampNs = timestampNs;
			btree.put(assignments, info);
			lowerIndex.add(lowerEnergy, assignments);
			upperIndex.add(upperEnergy, assignments);
		}

		public void setLowerBound(int[] assignments, double energy, long timestampNs) {
			ConfInfo info = btree.get(assignments);
			if (info == null) {
				info = new ConfInfo();
			} else {
				// remove old energy index entries if needed
				if (info.lowerTimestampNs != 0L) {
					lowerIndex.remove(info.lowerEnergy, assignments);
				}
			}
			info.lowerEnergy = energy;
			info.lowerTimestampNs = timestampNs;
			btree.put(assignments, info);
			lowerIndex.add(energy, assignments);
		}

		public void setUpperBound(int[] assignments, double energy, long timestampNs) {
			ConfInfo info = btree.get(assignments);
			if (info == null) {
				info = new ConfInfo();
			} else {
				// remove old energy index entries if needed
				if (info.upperTimestampNs != 0L) {
					upperIndex.remove(info.upperEnergy, assignments);
				}
			}
			info.upperEnergy = energy;
			info.upperTimestampNs = timestampNs;
			btree.put(assignments, info);
			upperIndex.add(energy, assignments);
		}

		public Conf get(int[] assignments) {

			ConfInfo info = btree.get(assignments);
			if (info == null) {
				return null;
			}

			return new Conf(assignments, info);
		}

		public ConfSearch.ScoredConf getScored(int[] assignments) {

			ConfInfo info = btree.get(assignments);
			if (info == null) {
				return null;
			}

			return new ConfSearch.ScoredConf(
				assignments,
				info.lowerTimestampNs == 0L ? Double.NaN : info.lowerEnergy
			);
		}

		public ConfSearch.EnergiedConf getEnergied(ConfSearch.ScoredConf conf) {

			ConfInfo info = btree.get(conf.getAssignments());
			if (info == null || info.upperTimestampNs == 0L) {
				return null;
			}

			return new ConfSearch.EnergiedConf(conf, info.upperEnergy);
		}

		public ConfSearch.EnergiedConf getEnergied(int[] assignments) {

			ConfInfo info = btree.get(assignments);
			if (info == null) {
				return null;
			}

			return new ConfSearch.EnergiedConf(
				assignments,
				info.lowerTimestampNs == 0L ? Double.NaN : info.lowerEnergy,
				info.upperTimestampNs == 0L ? Double.NaN : info.upperEnergy
			);
		}

		public void remove(int[] assignments) {
			ConfInfo info = btree.get(assignments);
			if (info != null) {
				if (info.lowerTimestampNs != 0L) {
					lowerIndex.remove(info.lowerEnergy, assignments);
				}
				if (info.upperTimestampNs != 0L) {
					upperIndex.remove(info.upperEnergy, assignments);
				}
				btree.remove(assignments);
			}
		}

		@Override
		public Iterator<Conf> iterator() {
			return Streams.of(btree.entryIterator())
				.map((entry) -> new Conf(
						entry.getKey(),
						entry.getValue()
					)
				)
				.iterator();
		}

		public Iterable<ConfSearch.ScoredConf> scoredConfs(SortOrder sort) {
			switch (sort) {

				case Assignment:
					return () -> Streams.of(iterator())
						.map((conf) -> conf.toScoredConf())
						.iterator();

				case Score:
					return () -> Streams.of(lowerIndex.iterator())
						.map((entry) -> new ConfSearch.ScoredConf(entry.getValue(), entry.getKey()))
						.iterator();

				case Energy:
					return () -> Streams.of(upperIndex.iterator())
						.map((entry) -> getScored(entry.getValue()))
						.filter((conf) -> conf != null)
						.iterator();

				default:
					throw new UnpossibleError();
			}
		}

		public Iterable<ConfSearch.EnergiedConf> energiedConfs(SortOrder sort) {
			switch (sort) {

				case Assignment:
					return () -> Streams.of(iterator())
						.map((conf) -> conf.toEnergiedConf())
						.iterator();

				case Score:
					return () -> Streams.of(lowerIndex.iterator())
						.map((entry) -> getEnergied(entry.getValue()))
						.filter((conf) -> conf != null)
						.iterator();

				case Energy:
					return () -> Streams.of(upperIndex.iterator())
						.map((entry) -> getEnergied(entry.getValue()))
						.filter((conf) -> conf != null)
						.iterator();

				default:
					throw new UnpossibleError();
			}
		}

		public Iterable<Double> lowerBounds() {
			return () -> lowerIndex.btree.keyIterator();
		}

		public Iterable<Double> upperBounds() {
			return () -> upperIndex.btree.keyIterator();
		}

		public List<Conf> getConfsByLowerBound(double energy) {
			List<int[]> multiAssignments = lowerIndex.get(energy);
			if (multiAssignments == null) {
				return null;
			}
			return Streams.of(multiAssignments)
				.map((assignments) -> get(assignments))
				.collect(Collectors.toList());
		}

		public List<Conf> getConfsByUpperBound(double energy) {
			List<int[]> multiAssignments = upperIndex.get(energy);
			if (multiAssignments == null) {
				return null;
			}
			return Streams.of(multiAssignments)
				.map((assignments) -> get(assignments))
				.collect(Collectors.toList());
		}

		public long size() {
			return btree.sizeLong();
		}

		public void flush() {
			ConfDB.this.flush();
		}
	}

	public class SequenceDB extends ConfTable {

		public final Sequence sequence;

		public SequenceDB(Sequence sequence) {
			super(getSequenceId(sequence));
			this.sequence = sequence;
		}

		private SequenceInfo getInfo() {
			return sequences.get(sequence);
		}

		private void setInfo(SequenceInfo info) {
			sequences.put(sequence, info);
		}

		public double getLowerEnergyOfUnsampledConfs() {
			return getInfo().lowerEnergyOfUnsampledConfs;
		}

		public void setLowerEnergyOfUnsampledConfs(double val) {
			SequenceInfo info = getInfo();
			info.lowerEnergyOfUnsampledConfs = val;
			setInfo(info);
		}

		public void updateLowerEnergyOfUnsampledConfs(double val) {
			SequenceInfo info = getInfo();
			if (Double.isNaN(info.lowerEnergyOfUnsampledConfs) || val > info.lowerEnergyOfUnsampledConfs) {
				info.lowerEnergyOfUnsampledConfs = val;
			}
			setInfo(info);
		}
	}

	private class EnergyIndex implements Iterable<Map.Entry<Double,int[]>> {

		public final BTreeMap<Double,List<int[]>> btree;

		public EnergyIndex(String id) {
			this.btree = db.treeMap(id)
				.keySerializer(Serializer.DOUBLE)
				.valueSerializer(new MultiAssignmentsSerializer())
				.createOrOpen();
		}

		public List<int[]> get(double energy) {
			return btree.get(energy);
		}

		public void add(double energy, int[] assignments) {
			List<int[]> multiAssignments = get(energy);
			if (multiAssignments == null) {
				multiAssignments = Arrays.asList(assignments);
			} else {

				// if the assignments are already there, don't update anything
				if (hasAssignments(multiAssignments, assignments)) {
					return;
				}

				multiAssignments.add(assignments);
			}
			btree.put(energy, multiAssignments);
		}

		public void remove(double energy, int[] assignments) {
			List<int[]> multiAssignments = get(energy);
			if (multiAssignments != null && removeAssignments(multiAssignments, assignments)) {
				if (multiAssignments.size() == 0) {
					btree.remove(energy);
				} else {
					btree.put(energy, multiAssignments);
				}
			}
		}

		private class EnergiedAssignments implements Map.Entry<Double,int[]> {

			public final double energy;
			public final int[] assignments;

			public EnergiedAssignments(double energy, int[] assignments) {
				this.energy = energy;
				this.assignments = assignments;
			}

			@Override
			public Double getKey() {
				return energy;
			}

			@Override
			public int[] getValue() {
				return assignments;
			}

			@Override
			public int[] setValue(int[] value) {
				throw new UnsupportedOperationException();
			}
		}

		@Override
		public Iterator<Map.Entry<Double,int[]>> iterator() {
			return Streams.of(btree.entryIterator())
				.flatMap((entry) ->
					entry.getValue().stream()
						.map((assignments) ->
							(Map.Entry<Double,int[]>)new EnergiedAssignments(entry.getKey(), assignments)
						)
				)
				.iterator();
		}

		// sigh... java arrays don't implement equals() etc
		private boolean hasAssignments(List<int[]> multiAssignments, int[] assignments) {
			for (int[] a : multiAssignments) {
				if (Arrays.equals(a, assignments)) {
					return true;
				}
			}
			return false;
		}
		private boolean removeAssignments(List<int[]> multiAssignments, int[] assignments) {
			Iterator<int[]> iter = multiAssignments.iterator();
			while (iter.hasNext()) {
				int[] a = iter.next();
				if (Arrays.equals(a, assignments)) {
					iter.remove();
					return true;
				}
			}
			return false;
		}
	}

	public final SimpleConfSpace confSpace;
	public final File file;

	private final DB db;
	private final HTreeMap<Sequence,SequenceInfo> sequences;
	private final Map<Sequence,SequenceDB> sequenceDBs;
	private final IntEncoding assignmentEncoding;

	public ConfDB(SimpleConfSpace confSpace) {
		this(confSpace, null);
	}

	public ConfDB(SimpleConfSpace confSpace, File file) {

		this.confSpace = confSpace;
		this.file = file;

		// determine conf encoding
		int maxAssignment = 0;
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			for (SimpleConfSpace.ResidueConf resConf : pos.resConfs) {
				maxAssignment = Math.max(maxAssignment, resConf.index);
			}
		}
		assignmentEncoding = IntEncoding.get(maxAssignment + 1); // +1 for the shift to move the unassigned value (-1) to non-negative

		// MapDB serializer for Sequence
		SimpleSerializer<Sequence> sequenceSerializer = new SimpleSerializer<Sequence>(SimpleSerializer.DynamicSize) {

			@Override
			public void serialize(@NotNull DataOutput2 out, @NotNull Sequence sequence)
			throws IOException {
				out.writeUTF(getSequenceId(sequence));
			}

			@Override
			public Sequence deserialize(@NotNull DataInput2 in, int available)
			throws IOException {
				return makeSequenceFromId(in.readUTF());
			}

			@Override
			public int compare(Sequence a, Sequence b) {

				// short circuit
				if (a == b) {
					return 0;
				}

				// lexicographical comparison
				for (SeqSpace.Position pos : confSpace.seqSpace.positions) {
					SeqSpace.ResType aResType = a.get(pos);
					SeqSpace.ResType bResType = b.get(pos);
					if (aResType != null && bResType != null) {
						// both not null, safe to compare
						int val = aResType.compareTo(bResType);
						if (val != 0) {
							return val;
						}
					} else if (aResType == null && bResType != null) {
						// a null, but not b, assume a < b
						return -1;
					} else if (aResType != null) {
						// b null, but not a, assume a > b
						return 1;
					}
					// both null, continue to next pos
				}

				return 0;
			}
		};

		// MapDB serialzier for SequenceInfo
		final int infoSize = Double.BYTES;
		SimpleSerializer<SequenceInfo> infoSerializer = new SimpleSerializer<SequenceInfo>(infoSize) {

			@Override
			public void serialize(@NotNull DataOutput2 out, @NotNull SequenceInfo info)
			throws IOException {
				out.writeDouble(info.lowerEnergyOfUnsampledConfs);
			}

			@Override
			public SequenceInfo deserialize(@NotNull DataInput2 in, int available)
			throws IOException {
				return new SequenceInfo(in.readDouble());
			}
		};

		// open the DB
		if (file != null) {
			db = DBMaker.fileDB(file)
				.transactionEnable() // turn on wite-ahead log, so the db survives JVM crashes
				.fileMmapEnableIfSupported() // use memory-mapped files if possible (can be much faster)
				.closeOnJvmShutdown()
				.make();
		} else {
			db = DBMaker.memoryDB()
				.make();
		}
		sequences = db.hashMap("sequences")
			.keySerializer(sequenceSerializer)
			.valueSerializer(infoSerializer)
			.createOrOpen();
		sequenceDBs = new HashMap<>();
	}

	private String getSequenceId(Sequence sequence) {
		return String.join(":", () ->
			sequence.seqSpace.positions.stream()
				.map((pos) -> (CharSequence)sequence.get(pos).name)
				.iterator()
		);
	}

	private Sequence makeSequenceFromId(String id) {
		Sequence sequence = confSpace.makeUnassignedSequence();
		String[] resTypes = id.split(":");
		for (SeqSpace.Position pos : confSpace.seqSpace.positions) {
			sequence.set(pos.resNum, resTypes[pos.index]);
		}
		return sequence;
	}

	public long getNumSequences() {
		// Java API means we're stuck with int-sized values here
		//return sequences.getSize();
		// but we could have lots of sequences, so we want a long
		return sequences.sizeLong();
	}

	public Iterable<Sequence> getSequences() {
		return (Set<Sequence>)sequences.keySet();
	}

	public SequenceDB getSequence(Sequence sequence) {

		// make sure the sequence spaces match
		if (sequence.seqSpace != confSpace.seqSpace) {
			throw new IllegalArgumentException("this sequence is from a different sequence space than the sequence space used by this db");
		}

		SequenceDB sdb = sequenceDBs.get(sequence);
		if (sdb == null) {
			sdb = new SequenceDB(sequence);
			sequenceDBs.put(sequence, sdb);
			if (!sequences.containsKey(sequence)) {
				sequences.put(sequence, new SequenceInfo());
			}
		}
		return sdb;
	}

	public void flush() {
		// In write-ahead mode, we don't actually have any transactions,
		// so there's nothing to commit in the traditional sense.
		// So in this case, "commit" flushes write caches to disk
		db.commit();
	}

	public void close() {
		flush();
		for (ConfTable sdb : sequenceDBs.values()) {
			sdb.btree.close();
		}
		sequenceDBs.clear();
		db.close();
	}

	@Override
	public void clean() {
		close();
	}
}
