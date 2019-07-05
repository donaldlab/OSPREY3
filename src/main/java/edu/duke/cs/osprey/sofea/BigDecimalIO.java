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


import edu.duke.cs.osprey.tools.MathTools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;


/**
 * serializes BigDecimal instances to as fixed-length byte arrays
 * based on the MathContext precision
 *
 * supports MathTools pseudo-values for BigDecimal like -/+inf, NaN
 */
public interface BigDecimalIO {

	void write(DataOutput out, BigDecimal d) throws IOException;
	BigDecimal read(DataInput in) throws IOException;


	public static enum Type {

		Null(null),
		NaN(MathTools.BigNaN),
		NegInf(MathTools.BigNegativeInfinity),
		PosInf(MathTools.BigPositiveInfinity),
		Zero(BigDecimal.ZERO),
		One(BigDecimal.ONE),
		Val(null);

		private final byte code = (byte)ordinal();
		private final BigDecimal value;

		Type(BigDecimal value) {
			this.value = value;
		}

		public static Type get(BigDecimal d) {
			if (d == null) {
				return Null;
			} else if (d == MathTools.BigNaN) {
				return NaN;
			} else if (d == MathTools.BigNegativeInfinity) {
				return NegInf;
			} else if (d == MathTools.BigPositiveInfinity) {
				return PosInf;
			} else if (d == BigDecimal.ZERO) {
				return Zero;
			} else if (d == BigDecimal.ONE) {
				return One;
			} else {
				return Val;
			}
		}

		public static Type get(byte code) {
			return values()[code];
		}
	}

	public static class Variable implements BigDecimalIO {

		@Override
		public void write(DataOutput out, BigDecimal d)
		throws IOException {

			Type type = Type.get(d);
			out.writeByte(type.code);

			if (type == Type.Val) {
				byte[] bytes = d.unscaledValue().toByteArray();
				out.writeInt(bytes.length);
				out.write(bytes);
				out.writeInt(d.scale());
			}
		}

		@Override
		public BigDecimal read(DataInput in)
			throws IOException {

			Type type = Type.get(in.readByte());

			if (type == Type.Val) {
				byte[] bytes = new byte[in.readInt()];
				in.readFully(bytes);
				return new BigDecimal(
					new BigInteger(bytes),
					in.readInt()
				);
			} else {
				return type.value;
			}
		}
	}

	public static class Fixed implements BigDecimalIO {

		public final int numBytes;

		private final int unscaledValBytes;
		private final int noncodeBytes;

		private static final int otherBytes = 1 + Integer.BYTES;

		public Fixed(int numBytes) {
			this.numBytes = numBytes;
			this.unscaledValBytes = numBytes - otherBytes;
			this.noncodeBytes = numBytes - 1;
		}

		public Fixed(MathContext mathContext) {
			this(calcUnscaledValBytes(mathContext.getPrecision()) + otherBytes);
		}

		/**
		 * determine the max number of bytes needed to store each BigDecimal instance
		 * based on the precision of the MathContext
		 */
		private static int calcUnscaledValBytes(int precision) {

			MathContext mc = new MathContext(precision, RoundingMode.HALF_UP);

			// start with our favorite one-digit prime number
			BigDecimal d = new BigDecimal("7.0");

			// multiply numbers together until the size stabilizes for several iterations
			int size = 0;
			int numSame = 0;
			while (numSame < 10) {

				d = d.multiply(d, mc);

				int newsize = d.unscaledValue().toByteArray().length;
				if (newsize > size) {
					size = newsize;
					numSame = 0;
				} else {
					numSame++;
				}
			}

			return size;
		}

		private void writePad(DataOutput out, int size)
		throws IOException {
			for (int i=0; i<size; i++) {
				out.writeByte(0);
			}
		}

		@Override
		public void write(DataOutput out, BigDecimal d)
		throws IOException {

			Type type = Type.get(d);
			out.writeByte(type.code);

			if (type == Type.Val) {
				byte[] bytes = d.unscaledValue().toByteArray();
				assert (bytes.length <= unscaledValBytes);
				writePad(out, unscaledValBytes - bytes.length);
				out.write(bytes);
				out.writeInt(d.scale());
			} else {
				writePad(out, noncodeBytes);
			}
		}

		private void readPad(DataInput in, int size)
		throws IOException {
			for (int i=0; i<size; i++) {
				in.readByte();
			}
		}

		@Override
		public BigDecimal read(DataInput in)
		throws IOException {

			Type type = Type.get(in.readByte());

			if (type == Type.Val) {
				byte[] bytes = new byte[unscaledValBytes];
				in.readFully(bytes);
				return new BigDecimal(
					new BigInteger(bytes),
					in.readInt()
				);
			} else {
				readPad(in, noncodeBytes);
				return type.value;
			}
		}
	}
}
