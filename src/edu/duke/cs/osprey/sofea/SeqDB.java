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
import edu.duke.cs.osprey.tools.*;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import org.jetbrains.annotations.NotNull;
import org.mapdb.*;
import org.mapdb.serializer.GroupSerializerObjectArray;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.*;
import java.util.function.Consumer;


public class SeqDB implements AutoCloseable {

	private static boolean isEmptySum(BigDecimalBounds z) {
		return MathTools.isZero(z.lower) && MathTools.isZero(z.upper);
	}

	private static BigDecimalBounds makeEmptySum() {
		return new BigDecimalBounds(BigDecimal.ZERO, BigDecimal.ZERO);
	}

	public static class SeqInfo {

		public final BigDecimalBounds[] zSumBounds;

		public SeqInfo(int size) {
			this.zSumBounds = new BigDecimalBounds[size];
		}

		public void setEmpty() {
			for (int i=0; i<zSumBounds.length; i++) {
				zSumBounds[i] = makeEmptySum();
			}
		}

		public boolean isEmpty() {
			for (BigDecimalBounds b : zSumBounds) {
				if (!isEmptySum(b)) {
					return false;
				}
			}
			return true;
		}

		public BigDecimalBounds get(MultiStateConfSpace.State state) {
			return zSumBounds[state.sequencedIndex];
		}

		@Override
		public int hashCode() {
			return HashCalculator.combineObjHashes(zSumBounds);
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof SeqInfo && equals((SeqInfo)other);
		}

		public boolean equals(SeqInfo other) {
			return Arrays.equals(this.zSumBounds, other.zSumBounds);
		}

		@Override
		public String toString() {
			return Streams.joinToString(zSumBounds, ", ", b -> Log.formatBigLn(b));
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
				s.serialize(out, data.zSumBounds[i]);
			}
		}

		@Override
		public SeqInfo deserialize(@NotNull DataInput2 in, int available)
		throws IOException {
			SeqInfo data = new SeqInfo(numBounds);
			for (int i=0; i<numBounds; i++) {
				data.zSumBounds[i] = s.deserialize(in, available);
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
	private final HTreeMap<int[],SeqInfo> sequencedSums;
	private final HTreeMap<Integer,BigDecimalBounds> unsequencedSums;

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
				.fileMmapEnableIfSupported() // use memory-mapped files if possible (can be much faster)
				// TODO: optimize allocation scheme? https://jankotek.gitbooks.io/mapdb/content/performance/
				.closeOnJvmShutdown()
				.make();
		} else {
			db = DBMaker.memoryDB()
				.make();
		}

		// open the tables

		sequencedSums = db.hashMap("sequenced-sums")
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

		unsequencedSums = db.hashMap("unsequenced-sums")
			.keySerializer(Serializer.INTEGER)
			.valueSerializer(new BigDecimalBoundsSerializer())
			.createOrOpen();
	}

	public BigMath bigMath() {
		return new BigMath(mathContext);
	}


	public class Transaction {

		private final Map<Sequence,SeqInfo> sequencedSums = new HashMap<>();
		private final Map<Integer,BigDecimalBounds> unsequencedSums = new HashMap<>();
		private boolean isEmpty = true;

		private Transaction() {
			// keep the constructor private
		}

		public BigMath bigMath() {
			return SeqDB.this.bigMath();
		}

		private void updateZSumBounds(MultiStateConfSpace.State state, Sequence seq, Consumer<BigDecimalBounds> f) {

			if (state.isSequenced) {

				// get the tx seq info, or empty sums
				SeqInfo seqInfo = sequencedSums.get(seq);
				if (seqInfo == null) {
					seqInfo = new SeqInfo(confSpace.sequencedStates.size());
					seqInfo.setEmpty();
				}

				f.accept(seqInfo.get(state));
				sequencedSums.put(seq, seqInfo);

			} else {

				// get the tx sum, or empty
				BigDecimalBounds sum = unsequencedSums.get(state.unsequencedIndex);
				if (sum == null) {
					sum = new BigDecimalBounds(BigDecimal.ZERO, BigDecimal.ZERO);
				}

				f.accept(sum);
				unsequencedSums.put(state.unsequencedIndex, sum);

			}

			isEmpty = false;
		}

		public void addZPath(MultiStateConfSpace.State state, Sequence seq, BigExp zPath, BigExp zSumUpper) {

			if (!zPath.isFinite() || !zSumUpper.isFinite()) {
				throw new IllegalArgumentException("Z must be finite: " + zPath + ", " + zSumUpper);
			}

			updateZSumBounds(state, seq, sum -> {
				sum.lower = bigMath()
					.set(sum.lower)
					.add(zPath.toBigDecimal(mathContext))
					.get();
				sum.upper = bigMath()
					.set(sum.upper)
					.add(zPath.toBigDecimal(mathContext))
					.sub(zSumUpper.toBigDecimal(mathContext))
					.get();
			});
		}

		public void addZSumUpper(MultiStateConfSpace.State state, Sequence seq, BigExp zSumUpper) {

			if (!zSumUpper.isFinite()) {
				throw new IllegalArgumentException("Z must be finite: " + zSumUpper);
			}

			updateZSumBounds(state, seq, sum ->
				sum.upper = bigMath()
					.set(sum.upper)
					.add(zSumUpper.toBigDecimal(mathContext))
					.get()
			);
		}

		public void subZSumUpper(MultiStateConfSpace.State state, Sequence seq, BigExp zSumUpper) {

			if (!zSumUpper.isFinite()) {
				throw new IllegalArgumentException("Z must be finite: " + zSumUpper);
			}

			updateZSumBounds(state, seq, sum ->
				sum.upper = bigMath()
					.set(sum.upper)
					.sub(zSumUpper.toBigDecimal(mathContext))
					.get()
			);
		}

		public boolean isEmpty() {
			return isEmpty;
		}

		private void combineSums(BigDecimalBounds sum, BigDecimalBounds oldSum) {
			sum.upper = bigMath()
				.set(sum.upper)
				.add(oldSum.upper)
				.get();
			sum.lower = bigMath()
				.set(sum.lower)
				.add(oldSum.lower)
				.get();
		}

		private void fixRoundoffError(BigDecimalBounds z) {
			// trust the lower bound more, since it's based on minimizations
			if (!z.isValid()) {
				// TODO: throw an Exception if the error is bigger than what we'd expect from roundoff?
				z.upper = z.lower;
			}
		}

		public void commit() {

			// short circuit
			if (isEmpty) {
				return;
			}

			// push writes to the db
			for (Map.Entry<Sequence,SeqInfo> entry : sequencedSums.entrySet()) {
				Sequence seq = entry.getKey();
				SeqInfo seqInfo = entry.getValue();

				// combine with the old sums if needed
				SeqInfo oldSeqInfo = SeqDB.this.sequencedSums.get(seq.rtIndices);
				if (oldSeqInfo != null) {
					for (MultiStateConfSpace.State state : confSpace.sequencedStates) {
						BigDecimalBounds sum = seqInfo.zSumBounds[state.sequencedIndex];
						BigDecimalBounds oldSum = oldSeqInfo.zSumBounds[state.sequencedIndex];
						combineSums(sum, oldSum);
						fixRoundoffError(sum);
					}
				}

				SeqDB.this.sequencedSums.put(seq.rtIndices, seqInfo);
			}

			for (Map.Entry<Integer,BigDecimalBounds> entry : unsequencedSums.entrySet()) {
				int unsequencedIndex = entry.getKey();
				BigDecimalBounds sum = entry.getValue();

				// combine with the old sum if needed
				BigDecimalBounds oldSum = SeqDB.this.unsequencedSums.get(unsequencedIndex);
				if (oldSum != null) {
					combineSums(sum, oldSum);
					fixRoundoffError(sum);
				}

				SeqDB.this.unsequencedSums.put(unsequencedIndex, sum);
			}

			db.commit();

			// reset state
			sequencedSums.clear();
			unsequencedSums.clear();
			isEmpty = true;
		}
	}

	public Transaction transaction() {
		return new Transaction();
	}

	@Override
	public void close() {
		db.close();
	}

	/**
	 * returns accumulated Z values for the queried state
	 * (you probably don't want this unless you're debugging)
	 */
	public BigDecimalBounds getUnsequencedSum(MultiStateConfSpace.State state) {
		BigDecimalBounds z = unsequencedSums.get(state.unsequencedIndex);
		if (z == null) {
			z = makeEmptySum();
		}
		return z;
	}

	/**
	 * returns the current Z bounds for the queried state
	 */
	public BigDecimalBounds getUnsequencedZSumBounds(MultiStateConfSpace.State state) {
		BigDecimalBounds z = unsequencedSums.get(state.unsequencedIndex);
		if (z == null) {
			z = new BigDecimalBounds(BigDecimal.ZERO, MathTools.BigPositiveInfinity);
		}
		return z;
	}

	/**
	 * returns accumulated state Z values for the queried sequence
	 * does not include Z uncertainty from ancestral partial sequences
	 * (you probably don't want this unless you're debugging)
	 */
	public SeqInfo getSequencedSums(Sequence seq) {
		SeqInfo seqInfo = sequencedSums.get(seq.rtIndices);
		if (seqInfo == null) {
			seqInfo = new SeqInfo(confSpace.sequencedStates.size());
			seqInfo.setEmpty();
		}
		return seqInfo;
	}

	/**
	 * returns accumulated state Z values for all sequences
	 */
	public Iterable<Map.Entry<Sequence,SeqInfo>> getSequencedSums() {

		return () -> new Iterator<Map.Entry<Sequence,SeqInfo>>() {

			Iterator<Map.Entry<int[],SeqInfo>> iter = sequencedSums.getEntries().iterator();

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
	 * returns the current state Z bounds for the queried sequence
	 * bounds for partial sequences only describe *unexplored* subtrees
	 * as more subtrees get explored, those Z values will be transfered to more fully-assigned sequences
	 * bounds for fully-explored partial sequences will be zero
	 */
	public SeqInfo getSequencedZSumBounds(Sequence seq) {
		SeqInfo seqInfo = getSequencedSums(seq);
		if (seq.isFullyAssigned()) {
			addZAncestry(seq, seqInfo);
		}
		return seqInfo;
	}

	/**
	 * returns Z bounds for all sequences
	 * returns bounds for both full and partial sequences
	 */
	public Iterable<Map.Entry<Sequence,SeqInfo>> getSequencedZSumBounds() {
		return () -> new Iterator<Map.Entry<Sequence,SeqInfo>>() {

			Iterator<Map.Entry<int[],SeqInfo>> iter = sequencedSums.getEntries().iterator();

			@Override
			public boolean hasNext() {
				return iter.hasNext();
			}

			@Override
			public Map.Entry<Sequence,SeqInfo> next() {

				Map.Entry<int[],SeqInfo> entry = iter.next();

				Sequence seq = new Sequence(confSpace.seqSpace, entry.getKey());
				SeqInfo seqInfo = entry.getValue();

				if (seq.isFullyAssigned()) {
					addZAncestry(seq, seqInfo);
				}

				return new AbstractMap.SimpleEntry<>(seq, seqInfo);
			}
		};
	}

	private void addZAncestry(Sequence seq, SeqInfo seqInfo) {

		int[] rtIndices = new int[confSpace.seqSpace.positions.size()];
		System.arraycopy(seq.rtIndices, 0, rtIndices, 0, rtIndices.length);

		// add uncertainty from partial sequence ancestry
		// NOTE: assumes tree pos order follows seq pos order
		for (int i=rtIndices.length - 1; i>=0; i--) {
			rtIndices[i] = Sequence.Unassigned;

			SeqInfo parentSeqInfo = sequencedSums.get(rtIndices);
			if (parentSeqInfo != null) {
				// couldn't that unexplored subtree contain no confs for this seq?
				// NOTE: don't add the lower bounds, the subtree need not necessarily contain confs for this sequence
				for (MultiStateConfSpace.State state : confSpace.sequencedStates) {
					seqInfo.zSumBounds[state.sequencedIndex].upper = bigMath()
						.set(seqInfo.zSumBounds[state.sequencedIndex].upper)
						.add(parentSeqInfo.zSumBounds[state.sequencedIndex].upper)
						.get();
				}
			}
		}
	}

	public String dump() {
		StringBuilder buf = new StringBuilder();
		buf.append("Unsequenced");
		for (MultiStateConfSpace.State state : confSpace.unsequencedStates) {
			BigDecimalBounds zSumBounds = getUnsequencedZSumBounds(state);
			buf.append(String.format("\n%10s  zSumBounds=%s w=%s",
				state.name,
				Log.formatBigLn(zSumBounds),
				Log.formatBigLn(zSumBounds.size(mathContext))
			));
		}
		buf.append("\nSequenced");
		for (Map.Entry<Sequence,SeqInfo> entry : getSequencedZSumBounds()) {
			Sequence seq = entry.getKey();
			SeqInfo seqInfo = entry.getValue();
			buf.append(String.format("\n[%s]", seq));
			for (MultiStateConfSpace.State state : confSpace.sequencedStates) {
				buf.append(String.format("\n%10s  zSumBounds=%s w=%s",
					state.name,
					Log.formatBigLn(seqInfo.get(state)),
					Log.formatBigLn(seqInfo.get(state).size(mathContext))
				));
			}
		}
		return buf.toString();
	}
}
