package edu.duke.cs.osprey.sofea;

import static edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.tools.*;
import org.jetbrains.annotations.NotNull;
import org.mapdb.*;
import org.mapdb.serializer.GroupSerializerObjectArray;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.AbstractMap;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.function.Consumer;


public class SeqDB implements AutoCloseable {

	public static class SeqInfo {

		public final BigDecimalBounds[] bounds;

		public SeqInfo(int numBounds) {
			this.bounds = new BigDecimalBounds[numBounds];
		}

		public void setZero() {
			for (int i=0; i<bounds.length; i++) {
				bounds[i] = new BigDecimalBounds(BigDecimal.ZERO, BigDecimal.ZERO);
			}
		}

		public boolean isZero() {
			for (BigDecimalBounds b : bounds) {
				if (!MathTools.isZero(b.lower) || !MathTools.isZero(b.upper)) {
					return false;
				}
			}
			return true;
		}

		@Override
		public int hashCode() {
			return HashCalculator.combineObjHashes(bounds);
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof SeqInfo && equals((SeqInfo)other);
		}

		public boolean equals(SeqInfo other) {
			return Arrays.equals(this.bounds, other.bounds);
		}
	}

	private static abstract class SimpleSerializer<T> extends GroupSerializerObjectArray<T> {

		public final int fixedSize;

		protected SimpleSerializer() {
			// dynamic size, rather than fixed
			this(-1);
		}

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
	}

	// NOTE: all int values get +1 when serialized, so we can accomodate the range [-1,maxVal]
	private static class IntArraySerializer extends SimpleSerializer<int[]> {

		private final IntEncoding encoding;
		private final int numPos;

		public IntArraySerializer(int maxVal, int numPos) {
			this(IntEncoding.get(maxVal + 1), numPos);
		}

		private IntArraySerializer(IntEncoding encoding, int numPos) {
			super(encoding.numBytes*numPos);
			this.encoding = encoding;
			this.numPos = numPos;
		}

		@Override
		public void serialize(@NotNull DataOutput2 out, @NotNull int[] data)
		throws IOException {
			assert (data.length == numPos);
			for (int i=0; i<numPos; i++) {
				encoding.write(out, data[i] + 1);
			}
		}

		@Override
		public int[] deserialize(@NotNull DataInput2 in, int available)
		throws IOException {
			int[] data = new int[numPos];
			for (int i=0; i<numPos; i++) {
				data[i] = encoding.read(in) - 1;
			}
			return data;
		}

		@Override
		public int compare(int[] a, int[] b) {
			throw new UnsupportedOperationException();
		}

		@Override
		public boolean equals(int[] a, int[] b) {
			return Arrays.equals(a, b);
		}

		@Override
		public int hashCode(@NotNull int[] data, int seed) {
			return DataIO.intHash(Arrays.hashCode(data) + seed);
		}
	}

	private static class BigDecimalSerializer extends SimpleSerializer<BigDecimal> {

		private BigDecimalIO io = new BigDecimalIO.Variable();

		@Override
		public void serialize(@NotNull DataOutput2 out, @NotNull BigDecimal data)
		throws IOException {
			io.write(out, data);
		}

		@Override
		public BigDecimal deserialize(@NotNull DataInput2 in, int available)
		throws IOException {
			return io.read(in);
		}

		@Override
		public int compare(BigDecimal a, BigDecimal b) {
			return MathTools.compare(a, b);
		}

		@Override
		public boolean equals(BigDecimal a, BigDecimal b) {
			return MathTools.isSameValue(a, b);
		}

		@Override
		public int hashCode(@NotNull BigDecimal data, int seed) {
			return DataIO.intHash(data.hashCode() + seed);
		}
	}

	private static class BigDecimalBoundsSerializer extends SimpleSerializer<BigDecimalBounds> {

		private final BigDecimalSerializer s = new BigDecimalSerializer();

		@Override
		public void serialize(@NotNull DataOutput2 out, @NotNull BigDecimalBounds data)
		throws IOException {
			s.serialize(out, data.lower);
			s.serialize(out, data.upper);
		}

		@Override
		public BigDecimalBounds deserialize(@NotNull DataInput2 in, int available)
		throws IOException {
			return new BigDecimalBounds(
				s.deserialize(in, available),
				s.deserialize(in, available)
			);
		}

		@Override
		public int compare(BigDecimalBounds a, BigDecimalBounds b) {
			throw new UnsupportedOperationException();
		}

		@Override
		public boolean equals(BigDecimalBounds a, BigDecimalBounds b) {
			return a.equals(b);
		}

		@Override
		public int hashCode(@NotNull BigDecimalBounds data, int seed) {
			return DataIO.intHash(data.hashCode() + seed);
		}
	}

	private static class SeqInfoSerializer extends SimpleSerializer<SeqInfo> {

		private final int numBounds;

		private final BigDecimalBoundsSerializer s = new BigDecimalBoundsSerializer();

		public SeqInfoSerializer(int numBounds) {
			this.numBounds = numBounds;
		}

		@Override
		public void serialize(@NotNull DataOutput2 out, @NotNull SeqInfo data)
		throws IOException {
			for (int i=0; i<numBounds; i++) {
				s.serialize(out, data.bounds[i]);
			}
		}

		@Override
		public SeqInfo deserialize(@NotNull DataInput2 in, int available)
		throws IOException {
			SeqInfo data = new SeqInfo(numBounds);
			for (int i=0; i<numBounds; i++) {
				data.bounds[i] = s.deserialize(in, available);
			}
			return data;
		}

		@Override
		public int compare(SeqInfo a, SeqInfo b) {
			throw new UnsupportedOperationException();
		}

		@Override
		public boolean equals(SeqInfo a, SeqInfo b) {
			return a.equals(b);
		}

		@Override
		public int hashCode(@NotNull SeqInfo data, int seed) {
			return DataIO.intHash(data.hashCode() + seed);
		}
	}


	public final MultiStateConfSpace confSpace;
	public final MathContext mathContext;
	public final File file;

	private final DB db;
	private final HTreeMap<int[],SeqInfo> sequencedBounds;
	private final HTreeMap<Integer,BigDecimalBounds> unsequencedBounds;

	public SeqDB(MultiStateConfSpace confSpace, MathContext mathContext) {
		this(confSpace, mathContext, null);
	}

	public SeqDB(MultiStateConfSpace confSpace, MathContext mathContext, File file) {

		this.confSpace = confSpace;
		this.mathContext = mathContext;
		this.file = file;

		// open the DB
		if (file != null) {
			db = DBMaker.fileDB(file)
				.transactionEnable()
				.fileMmapEnableIfSupported() // use memory-mapped files if possible (can be much faster)
				// TODO: optimize allocation scheme? https://jankotek.gitbooks.io/mapdb/content/performance/
				.closeOnJvmShutdown()
				.make();
		} else {
			db = DBMaker.memoryDB()
				.transactionEnable()
				.make();
		}

		// open the tables

		sequencedBounds = db.hashMap("sequenced-bounds")
			.keySerializer(new IntArraySerializer(
				confSpace.seqSpace.positions.stream()
					.flatMap(pos -> pos.resTypes.stream())
					.mapToInt(resType -> resType.index)
					.max()
					.orElse(-1),
				confSpace.seqSpace.positions.size()
			))
			.valueSerializer(new SeqInfoSerializer(confSpace.sequencedStates.size()))
			.createOrOpen();

		unsequencedBounds = db.hashMap("unsequenced-bounds")
			.keySerializer(Serializer.INTEGER)
			.valueSerializer(new BigDecimalBoundsSerializer())
			.createOrOpen();
	}

	private BigMath bigMath() {
		return new BigMath(mathContext);
	}

	private Sequence makeSeq(MultiStateConfSpace.State state, ConfIndex index) {
		Sequence seq = confSpace.seqSpace.makeUnassignedSequence();
		for (int i=0; i<index.numDefined; i++) {
			SimpleConfSpace.Position confPos = state.confSpace.positions.get(index.definedPos[i]);
			SimpleConfSpace.ResidueConf resConf = confPos.resConfs.get(index.definedRCs[i]);
			if (confPos.seqPos != null) {
				seq.set(confPos.resNum, resConf.template.name);
			}
		}
		return seq;
	}

	private void updateSeqZ(MultiStateConfSpace.State state, ConfIndex index, Consumer<BigDecimalBounds> f) {

		if (state.isSequenced) {

			Sequence seq = makeSeq(state, index);

			// get the existing bounds, or a [0,0] bound
			SeqInfo seqInfo = sequencedBounds.get(seq.rtIndices);
			if (seqInfo == null) {
				seqInfo = new SeqInfo(confSpace.sequencedStates.size());
				seqInfo.setZero();
			}

			// apply the modification to the state bounds
			BigDecimalBounds bounds = seqInfo.bounds[state.sequencedIndex];
			f.accept(bounds);

			// did we end up with [0,0] again?
			if (seqInfo.isZero()) {

				// yup, just delete the whole record. no need to store [0,0] bounds
				sequencedBounds.remove(seq.rtIndices);

			} else {

				// otherwise, update the db
				sequencedBounds.put(seq.rtIndices, seqInfo);
			}

		} else {

			BigDecimalBounds bounds = unsequencedBounds.get(state.index);
			if (bounds == null) {
				bounds = new BigDecimalBounds(BigDecimal.ZERO, BigDecimal.ZERO);
			}

			// apply the modification
			f.accept(bounds);

			// did we end up with [0,0] again?
			if (MathTools.isZero(bounds.lower) && MathTools.isZero(bounds.upper)) {

				// yup, just delete the whole record. no need to store [0,0] bounds
				unsequencedBounds.remove(state.index);

			} else {

				// otherwise, update the db
				unsequencedBounds.put(state.index, bounds);
			}
		}
	}

	public void addSeqZ(MultiStateConfSpace.State state, ConfIndex index, BigDecimal Z) {
		updateSeqZ(state, index, bounds -> {
			bounds.lower = bigMath()
				.set(bounds.lower)
				.add(Z)
				.get();
			bounds.upper = bigMath()
				.set(bounds.upper)
				.add(Z)
				.get();
		});
	}

	public void addSeqZ(MultiStateConfSpace.State state, ConfIndex index, BigDecimalBounds Z) {
		updateSeqZ(state, index, bounds -> {
			bounds.lower = bigMath()
				.set(bounds.lower)
				.add(Z.lower)
				.get();
			bounds.upper = bigMath()
				.set(bounds.upper)
				.add(Z.upper)
				.get();
		});
	}

	public void subSeqZ(MultiStateConfSpace.State state, ConfIndex index, BigDecimalBounds Z) {
		updateSeqZ(state, index, bounds -> {
			bounds.upper = bigMath()
				.set(bounds.upper)
				.sub(Z.upper)
				.atLeast(0.0) // NOTE: roundoff error can cause this to drop below 0
				.get();
			bounds.lower = bigMath()
				.set(bounds.lower)
				.sub(Z.lower)
				.atMost(bounds.upper) // don't exceed the upper bound due to roundoff error
				.get();
		});
	}

	public void commit() {
		db.commit();
	}

	public void rollback() {
		db.rollback();
	}

	@Override
	public void close() {
		db.close();
	}

	public BigDecimalBounds getUnsequenced(int stateIndex) {
		return unsequencedBounds.get(stateIndex);
	}

	/**
	 * returns accumulated state Z values for all sequences
	 * does not include Z uncertainty from ancestral partial sequences
	 * (you probably don't want this unless you're debugging)
	 */
	public Iterable<Map.Entry<Sequence,SeqInfo>> getSequencesWithoutAncestry() {

		return () -> new Iterator<Map.Entry<Sequence,SeqInfo>>() {

			Iterator<Map.Entry<int[],SeqInfo>> iter = sequencedBounds.getEntries().iterator();

			@Override
			public boolean hasNext() {
				return iter.hasNext();
			}

			@Override
			public Map.Entry<Sequence, SeqInfo> next() {

				Map.Entry<int[],SeqInfo> entry = iter.next();

				Sequence seq = new Sequence(confSpace.seqSpace, entry.getKey());
				SeqInfo seqInfo = entry.getValue();

				return new AbstractMap.SimpleEntry<>(seq, seqInfo);
			}
		};
	}

	/**
	 * returns the current state Z bounds for all sequences
	 */
	public Iterable<Map.Entry<Sequence,SeqInfo>> getSequences() {

		int[] rtIndices = new int[confSpace.seqSpace.positions.size()];

		return () -> new Iterator<Map.Entry<Sequence,SeqInfo>>() {

			Iterator<Map.Entry<int[],SeqInfo>> iter = sequencedBounds.getEntries().iterator();

			@Override
			public boolean hasNext() {
				return iter.hasNext();
			}

			@Override
			public Map.Entry<Sequence, SeqInfo> next() {

				Map.Entry<int[],SeqInfo> entry = iter.next();

				Sequence seq = new Sequence(confSpace.seqSpace, entry.getKey());
				SeqInfo seqInfo = entry.getValue();

				if (seq.isFullyAssigned()) {

					// add uncertainty from partial sequence ancestry
					// NOTE: assumes tree pos order follows seq pos order
					System.arraycopy(seq.rtIndices, 0, rtIndices, 0, rtIndices.length);
					for (int i=rtIndices.length - 1; i>=0; i--) {
						rtIndices[i] = Sequence.Unassigned;

						SeqInfo parentSeqInfo = sequencedBounds.get(rtIndices);
						if (parentSeqInfo != null) {
							// couldn't that unexplored subtree contain no confs for this seq?
							// NOTE: don't add the lower bounds, the subtree need not necessarily contain confs for this sequence
							for (MultiStateConfSpace.State state : confSpace.sequencedStates) {
								seqInfo.bounds[state.sequencedIndex].upper = bigMath()
									.set(seqInfo.bounds[state.sequencedIndex].upper)
									.add(parentSeqInfo.bounds[state.sequencedIndex].upper)
									.get();
							}
						}
					}
				}

				return new AbstractMap.SimpleEntry<>(seq, seqInfo);
			}
		};
	}
}
