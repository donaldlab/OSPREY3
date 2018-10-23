package edu.duke.cs.osprey.sofea;

import static edu.duke.cs.osprey.tools.Log.log;
import static edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.IntEncoding;
import edu.duke.cs.osprey.tools.MathTools;
import org.jetbrains.annotations.NotNull;
import org.mapdb.*;
import org.mapdb.serializer.GroupSerializerObjectArray;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.Arrays;
import java.util.Map;
import java.util.function.Consumer;


public class SeqDB implements AutoCloseable {

	public class SeqInfo {

		public final Sequence seq;
		public final BigDecimalBounds bounds;

		public SeqInfo(Sequence seq, BigDecimalBounds bounds) {
			this.seq = seq;
			this.bounds = bounds;
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

		@Override
		public void serialize(@NotNull DataOutput2 out, @NotNull BigDecimal data)
		throws IOException {
			if (data == MathTools.BigNaN) {
				out.writeByte(1);
			} else if (data == MathTools.BigNegativeInfinity) {
				out.writeByte(2);
			} else if (data == MathTools.BigPositiveInfinity) {
				out.writeByte(3);
			} else if (data == BigDecimal.ZERO) {
				out.writeByte(4);
			} else if (data == BigDecimal.ONE) {
				out.writeByte(5);
			} else {
				out.writeByte(6);
				byte[] bytes = data.unscaledValue().toByteArray();
				// TODO: could use short or byte here?
				out.writeInt(bytes.length);
				out.write(bytes);
				out.writeInt(data.scale());
			}
		}

		@Override
		public BigDecimal deserialize(@NotNull DataInput2 in, int available)
		throws IOException {
			byte type = in.readByte();
			switch (type) {
				case 1: return MathTools.BigNaN;
				case 2: return MathTools.BigNegativeInfinity;
				case 3: return MathTools.BigPositiveInfinity;
				case 4: return BigDecimal.ZERO;
				case 5: return BigDecimal.ONE;
				case 6: {
					byte[] bytes = new byte[in.readInt()];
					in.readFully(bytes);
					return new BigDecimal(
						new BigInteger(bytes),
						in.readInt()
					);
				}
				default: throw new IOException("unrecognized BigDecimal type: " + type);
			}
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


	public final SimpleConfSpace confSpace;
	public final MathContext mathContext;
	public final File file;

	private final DB db;
	private final HTreeMap<int[],BigDecimalBounds> sequences;

	public SeqDB(SimpleConfSpace confSpace, MathContext mathContext) {
		this(confSpace, mathContext, null);
	}

	public SeqDB(SimpleConfSpace confSpace, MathContext mathContext, File file) {

		this.confSpace = confSpace;
		this.mathContext = mathContext;
		this.file = file;

		// open the DB
		if (file != null) {
			db = DBMaker.fileDB(file)
				.transactionEnable() // turn on wite-ahead log, so the db survives JVM crashes
				.fileMmapEnableIfSupported() // use memory-mapped files if possible (can be much faster)
				// TODO: optimize allocation scheme? https://jankotek.gitbooks.io/mapdb/content/performance/
				.closeOnJvmShutdown()
				.make();
		} else {
			db = DBMaker.memoryDB()
				.make();
		}

		// open the tables
		sequences = db.hashMap("sequences")
			.keySerializer(new IntArraySerializer(
				confSpace.seqSpace.positions.stream()
					.flatMap(pos -> pos.resTypes.stream())
					.mapToInt(resType -> resType.index)
					.max()
					.orElse(-1),
				confSpace.positions.size()
			))
			.valueSerializer(new BigDecimalBoundsSerializer())
			.createOrOpen();
	}

	private BigMath bigMath() {
		return new BigMath(mathContext);
	}

	// TODO: aggregate writes in memory before hitting the DB?

	private Sequence makeSeq(ConfIndex index) {
		Sequence seq = confSpace.seqSpace.makeUnassignedSequence();
		for (int i=0; i<index.numDefined; i++) {
			SimpleConfSpace.Position confPos = confSpace.positions.get(index.definedPos[i]);
			SimpleConfSpace.ResidueConf resConf = confPos.resConfs.get(index.definedRCs[i]);
			if (confPos.seqPos != null) {
				seq.set(confPos.seqPos, resConf.template.name);
			}
		}
		return seq;
	}

	private void updateSeqZ(Sequence seq, Consumer<BigDecimalBounds> f) {

		// get the existing bounds, or a [0,0] bound
		BigDecimalBounds bounds = sequences.get(seq.rtIndices);
		if (bounds == null) {
			bounds = new BigDecimalBounds(BigDecimal.ZERO, BigDecimal.ZERO);
		}

		// apply the modification
		f.accept(bounds);

		// did we end up with [0,0] again?
		if (MathTools.isZero(bounds.lower) && MathTools.isZero(bounds.upper)) {

			// yup, just delete the whole record. no need to store [0,0] bounds
			sequences.remove(seq.rtIndices);

		} else {

			// otherwise, update the db
			sequences.put(seq.rtIndices, bounds);
		}
	}

	public void addSeqZ(ConfIndex index, BigDecimal Z) {
		updateSeqZ(makeSeq(index), bounds -> {
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

	public void addSeqZ(ConfIndex index, BigDecimalBounds Z) {
		updateSeqZ(makeSeq(index), bounds -> {
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

	public void subSeqZ(ConfIndex index, BigDecimalBounds Z) {
		updateSeqZ(makeSeq(index), bounds -> {
			bounds.lower = bigMath()
				.set(bounds.lower)
				.sub(Z.lower)
				.atLeast(0.0) // NOTE: roundoff error can cause this to drop below 0
				.get();
			bounds.upper = bigMath()
				.set(bounds.upper)
				.sub(Z.upper)
				.atLeast(0.0)
				.get();
		});
	}

	public void commit() {
		db.commit();
	}

	@Override
	public void close() {
		db.close();
	}

	// TEMP
	public void dump() {
		log("Seq DB raw:");
		for (Map.Entry<int[],BigDecimalBounds> entry : sequences.getEntries()) {
			Sequence seq = new Sequence(confSpace.seqSpace, entry.getKey());
			BigDecimalBounds bounds = entry.getValue();
			log("\t%s  [%s]", SofeaLab.dump(bounds), seq);
		}
	}

	public void dumpPartialSequences() {

		log("Seq DB: partial sequences, unexplored subtrees:");

		for (Map.Entry<int[],BigDecimalBounds> entry : sequences.getEntries()) {
			Sequence seq = new Sequence(confSpace.seqSpace, entry.getKey());

			if (seq.isFullyAssigned()) {
				continue;
			}

			BigDecimalBounds bounds = entry.getValue();
			log("\t%s  [%s]", SofeaLab.dump(bounds), seq);
		}
	}

	public void dumpSequences() {

		log("Seq DB: sequences:");

		int[] rtIndices = new int[confSpace.seqSpace.positions.size()];

		for (Map.Entry<int[],BigDecimalBounds> entry : sequences.getEntries()) {
			Sequence seq = new Sequence(confSpace.seqSpace, entry.getKey());

			if (!seq.isFullyAssigned()) {
				continue;
			}

			BigDecimalBounds bounds = entry.getValue();
			BigMath lomath = bigMath().set(bounds.lower);
			BigMath himath = bigMath().set(bounds.upper);

			// add uncertainty from partial sequence ancestry
			// NOTE: assumes tree pos order follows seq pos order
			System.arraycopy(seq.rtIndices, 0, rtIndices, 0, rtIndices.length);
			for (int i=rtIndices.length - 1; i>=0; i--) {
				rtIndices[i] = Sequence.Unassigned;

				BigDecimalBounds parentBounds = sequences.get(rtIndices);
				if (parentBounds != null) {
					lomath.add(parentBounds.lower);
					himath.add(parentBounds.upper);
				}
			}

			bounds = new BigDecimalBounds(lomath.get(), himath.get());
			log("\t%s  %.6f  [%s]", SofeaLab.dump(bounds), SofeaLab.getBoundsDelta(bounds), seq);
		}
	}
}
