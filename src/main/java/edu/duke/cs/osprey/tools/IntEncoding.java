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

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

public enum IntEncoding {

	// NOTE: this order is important for get()
	Byte(1, 255) {

		@Override
		public void write(DataOutput out, int val)
			throws IOException {
			out.writeByte(val);
		}

		@Override
		public int read(DataInput in)
			throws IOException {
			return in.readUnsignedByte();
		}
	},
	Short(2, 32767) {

		@Override
		public void write(DataOutput out, int val)
			throws IOException {
			out.writeShort(val);
		}

		@Override
		public int read(DataInput in)
			throws IOException {
			return in.readUnsignedShort();
		}
	},
	Int(4, Integer.MAX_VALUE) {

		@Override
		public void write(DataOutput out, int val)
			throws IOException {
			out.writeInt(val);
		}

		@Override
		public int read(DataInput in)
			throws IOException {
			return in.readInt();
		}
	};

	public final int numBytes;
	public final int maxValue;

	IntEncoding(int numBytes, int maxValue) {
		this.numBytes = numBytes;
		this.maxValue = maxValue;
	}

	public static IntEncoding get(int maxVal) {

		for (IntEncoding encoding : values()) {
			if (maxVal <= encoding.maxValue) {
				return encoding;
			}
		}

		// silly compiler... this can't happen, right?
		throw new Error("unpossible, unless I'm an idiot");
	}

	public abstract void write(DataOutput out, int val) throws IOException;
	public abstract int read(DataInput in) throws IOException;
}
