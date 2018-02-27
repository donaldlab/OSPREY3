package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.IntEncoding;
import edu.duke.cs.osprey.tools.Streams;
import org.jetbrains.annotations.NotNull;
import org.mapdb.*;
import org.mapdb.serializer.GroupSerializerObjectArray;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Consumer;
import java.util.function.Function;


public class ConfDB {

	public static interface UserWithReturn<T> {

		T use(ConfDB confdb) throws Exception;

		default T usePassExceptions(ConfDB confdb) {
			try {
				return use(confdb);
			} catch (Exception ex) {
				throw new RuntimeException(ex);
			}
		}
	}

	public static interface User {

		void use(ConfDB confdb) throws Exception;

		default void usePassExceptions(ConfDB confdb) {
			try {
				use(confdb);
			} catch (Exception ex) {
				throw new RuntimeException(ex);
			}
		}
	}

	public static <T> T useIfNeeded(SimpleConfSpace confSpace, File file, UserWithReturn<T> user) {

		if (confSpace == null || file == null) {

			// no DB? just pass null
			return user.usePassExceptions(null);

		} else {

			// open the db and make sure it gets cleaned up properly
			return new ConfDB(confSpace, file).use(user);
		}
	}

	public static void useIfNeeded(SimpleConfSpace confSpace, File file, User user) {

		if (confSpace == null || file == null) {

			// no DB? just pass null
			user.usePassExceptions(null);

		} else {

			// open the db and make sure it gets cleaned up properly
			new ConfDB(confSpace, file).use(user);
		}
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

		public ConfTable(String id) {

			int numPos = confSpace.positions.size();

			// MapDB serializer for assignments
			SimpleSerializer<int[]> assignmentsSerializer = new SimpleSerializer<int[]>(SimpleSerializer.DynamicSize) {

				@Override
				public void serialize(@NotNull DataOutput2 out, @NotNull int[] assignments)
				throws IOException {
					for (int i=0; i<numPos; i++) {
						assignmentEncoding.write(out, assignments[i]);
					}
				}

				@Override
				public int[] deserialize(@NotNull DataInput2 in, int available)
				throws IOException {
					int[] assignments = new int[numPos];
					for (int i=0; i<numPos; i++) {
						assignments[i] = assignmentEncoding.read(in);
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
			};

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
				.keySerializer(assignmentsSerializer)
				.valueSerializer(confInfoSerializer)
				.createOrOpen();
		}

		public void setBounds(int[] assignments, double lowerEnergy, double upperEnergy, long timestampNs) {
			btree.put(assignments, new ConfInfo(lowerEnergy, timestampNs, upperEnergy, timestampNs));
		}

		public void setLowerBound(int[] assignments, double energy, long timestampNs) {
			ConfInfo info = btree.get(assignments);
			if (info == null) {
				info = new ConfInfo();
			}
			info.lowerEnergy = energy;
			info.lowerTimestampNs = timestampNs;
			btree.put(assignments, info);
		}

		public void setUpperBound(int[] assignments, double energy, long timestampNs) {
			ConfInfo info = btree.get(assignments);
			if (info == null) {
				info = new ConfInfo();
			}
			info.upperEnergy = energy;
			info.upperTimestampNs = timestampNs;
			btree.put(assignments, info);
			db.commit();
		}

		public Conf get(int[] assignments) {

			ConfInfo info = btree.get(assignments);
			if (info == null) {
				return null;
			}

			return new Conf(assignments, info);
		}

		@NotNull
		@Override
		public Iterator<Conf> iterator() {
			Iterator<Map.Entry<int[],ConfInfo>> iter = btree.entryIterator();

			return new Iterator<Conf>() {

				@Override
				public boolean hasNext() {
					return iter.hasNext();
				}

				@Override
				public Conf next() {

					Map.Entry<int[],ConfInfo> entry = iter.next();

					return new Conf(
						entry.getKey(),
						entry.getValue()
					);
				}
			};
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

	public final SimpleConfSpace confSpace;
	public final File file;

	private final DB db;
	private final HTreeMap<Sequence,SequenceInfo> sequences;
	private final Map<Sequence,SequenceDB> sequenceDBs;
	private final IntEncoding assignmentEncoding;

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
		assignmentEncoding = IntEncoding.get(maxAssignment);

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
				for (SimpleConfSpace.Position pos : confSpace.positions) {
					String aResType = a.get(pos);
					String bResType = b.get(pos);
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
		db = DBMaker.fileDB(file)
			.transactionEnable() // turn on wite-ahead log, so the db survives JVM crashes
			.fileMmapEnableIfSupported() // use memory-mapped files if possible (can be much faster)
			.closeOnJvmShutdown()
			.make();
		sequences = db.hashMap("sequences")
			.keySerializer(sequenceSerializer)
			.valueSerializer(infoSerializer)
			.createOrOpen();
		sequenceDBs = new HashMap<>();
	}

	private String getSequenceId(Sequence sequence) {
		return String.join(":", Streams.toIterable(
			sequence.confSpace.positions.stream()
				.map((pos) -> sequence.get(pos))
		));
	}

	private Sequence makeSequenceFromId(String id) {
		Sequence sequence = Sequence.makeUnassigned(confSpace);
		String[] resTypes = id.split(":");
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			sequence.set(pos, resTypes[pos.index]);
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

		// make sure the conf spaces match
		if (sequence.confSpace != confSpace) {
			throw new IllegalArgumentException("this sequence is from a different conf space than the conf space used by this db");
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

	public <T> T use(UserWithReturn<T> user) {
		try {
			return user.usePassExceptions(this);
		} finally {
			close();
		}
	}

	public void use(User user) {
		try {
			user.usePassExceptions(this);
		} finally {
			close();
		}
	}
}
