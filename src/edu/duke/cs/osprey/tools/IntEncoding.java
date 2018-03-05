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
