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

package edu.duke.cs.osprey.tools;

import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.sofea.BigDecimalIO;
import org.jetbrains.annotations.NotNull;
import org.mapdb.DataIO;
import org.mapdb.DataInput2;
import org.mapdb.DataOutput2;
import org.mapdb.serializer.GroupSerializer;
import org.mapdb.serializer.GroupSerializerObjectArray;

import java.io.IOException;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class MapDBTools {

	public static abstract class SimpleSerializer<T> extends GroupSerializerObjectArray<T> {

		public static final int DynamicSize = -1;

		public final int fixedSize;

		protected SimpleSerializer() {
			// dynamic size, rather than fixed
			this(DynamicSize);
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
	public static class IntArraySerializer extends SimpleSerializer<int[]> {

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

	public static class ValuesSerializer<T> extends SimpleSerializer<List<T>> {

		private final GroupSerializer<T> valueSerializer;

		public ValuesSerializer(GroupSerializer<T> valueSerializer) {
			this.valueSerializer = valueSerializer;
		}

		@Override
		public void serialize(@NotNull DataOutput2 out, @NotNull List<T> data)
		throws IOException {
			out.writeInt(data.size());
			for (T value : data) {
				valueSerializer.serialize(out, value);
			}
		}

		@Override
		public List<T> deserialize(@NotNull DataInput2 in, int available)
		throws IOException {
			int size = in.readInt();
			var data = new ArrayList<T>(size);
			for (int i=0; i<size; i++) {
				data.add(valueSerializer.deserialize(in, available));
			}
			return data;
		}
	}

	public static class SequenceSerializer extends SimpleSerializer<Sequence> {

		public final SeqSpace seqSpace;

		public SequenceSerializer(SeqSpace seqSpace) {
			super(SimpleSerializer.DynamicSize);
			this.seqSpace = seqSpace;
		}

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
			for (SeqSpace.Position pos : seqSpace.positions) {
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

		public String getSequenceId(Sequence sequence) {
			return Streams.joinToString(
				sequence.seqSpace.positions,
				":",
				pos -> sequence.get(pos).name
			);
		}

		public Sequence makeSequenceFromId(String id) {
			Sequence sequence = seqSpace.makeUnassignedSequence();
			String[] resTypes = id.split(":");
			for (SeqSpace.Position pos : seqSpace.positions) {
				sequence.set(pos.resNum, resTypes[pos.index]);
			}
			return sequence;
		}
	}

	public static class BigDecimalSerializer extends SimpleSerializer<BigDecimal> {

		private final BigDecimalIO io = new BigDecimalIO.Variable();

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

	public static class BigDecimalBoundsSerializer extends SimpleSerializer<MathTools.BigDecimalBounds> {

		private final BigDecimalSerializer s = new BigDecimalSerializer();

		@Override
		public void serialize(@NotNull DataOutput2 out, @NotNull MathTools.BigDecimalBounds data)
		throws IOException {
			s.serialize(out, data.lower);
			s.serialize(out, data.upper);
		}

		@Override
		public MathTools.BigDecimalBounds deserialize(@NotNull DataInput2 in, int available)
		throws IOException {
			return new MathTools.BigDecimalBounds(
				s.deserialize(in, available),
				s.deserialize(in, available)
			);
		}

		@Override
		public int compare(MathTools.BigDecimalBounds a, MathTools.BigDecimalBounds b) {
			throw new UnsupportedOperationException();
		}

		@Override
		public boolean equals(MathTools.BigDecimalBounds a, MathTools.BigDecimalBounds b) {
			return a.equals(b);
		}

		@Override
		public int hashCode(@NotNull MathTools.BigDecimalBounds data, int seed) {
			return DataIO.intHash(data.hashCode() + seed);
		}
	}

	public static class BigExpSerializer extends MapDBTools.SimpleSerializer<BigExp> {

		public BigExpSerializer() {
			super(Double.BYTES + Integer.BYTES);
		}

		@Override
		public void serialize(@NotNull DataOutput2 out, @NotNull BigExp value)
		throws IOException {
			out.writeDouble(value.fp);
			out.writeInt(value.exp);
		}

		@Override
		public BigExp deserialize(@NotNull DataInput2 input, int available)
		throws IOException {
			return new BigExp(
				input.readDouble(),
				input.readInt()
			);
		}

		@Override
		public int compare(BigExp a, BigExp b) {
			return a.compareTo(b);
		}

		@Override
		public boolean equals(BigExp a, BigExp b) {
			return a.equals(b);
		}
	}
}
