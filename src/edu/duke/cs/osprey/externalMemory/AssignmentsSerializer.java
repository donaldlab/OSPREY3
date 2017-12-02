package edu.duke.cs.osprey.externalMemory;

import java.nio.ByteBuffer;
import java.util.Arrays;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.tpie.EntrySize;

public abstract class AssignmentsSerializer {

	public static class WrongEncodingException extends RuntimeException {

		private static final long serialVersionUID = -7973280483158957646L;

		public WrongEncodingException(Encoding encoding, int val) {
			super(String.format("%s encoding is wrong for value %d. This is definitely a bug",
				encoding.name(),
				val
			));
		}
	}

	public static class WrongNumberOfAssignmentsException extends RuntimeException {

		private static final long serialVersionUID = -2877773001328314147L;

		public WrongNumberOfAssignmentsException(RCs rcs, int[] assignments) {
			super(String.format("Serializer expected %d design positions, but got %d instead. (%s)",
				rcs.getNumPos(),
				assignments.length,
				Arrays.toString(assignments)
			));
		}
	}
	
	public static enum Encoding {
		
		// NOTE: this order is important for pickBest()
		
		Byte(java.lang.Byte.BYTES, java.lang.Byte.MAX_VALUE) {
			
			@Override
			public void write(int[] vals, ByteBuffer buf) {
				for (int i=0; i<vals.length; i++) {

					// just in case...
					if (vals[i] > java.lang.Byte.MAX_VALUE) {
						throw new WrongEncodingException(this, vals[i]);
					}

					buf.put((byte)vals[i]);
				}
			}
			
			@Override
			public void read(ByteBuffer buf, int[] vals) {
				for (int i=0; i<vals.length; i++) {
					vals[i] = buf.get();
				}
			}
		},
		Short(java.lang.Short.BYTES, java.lang.Short.MAX_VALUE) {
			
			@Override
			public void write(int[] vals, ByteBuffer buf) {
				for (int i=0; i<vals.length; i++) {

					// just in case...
					if (vals[i] > java.lang.Short.MAX_VALUE) {
						throw new WrongEncodingException(this, vals[i]);
					}

					buf.putShort((short)vals[i]);
				}
			}
			
			@Override
			public void read(ByteBuffer buf, int[] vals) {
				for (int i=0; i<vals.length; i++) {
					vals[i] = buf.getShort();
				}
			}
		},
		Int(java.lang.Integer.BYTES, java.lang.Integer.MAX_VALUE) {
			
			@Override
			public void write(int[] vals, ByteBuffer buf) {
				for (int i=0; i<vals.length; i++) {
					buf.putInt(vals[i]);
				}
			}
			
			@Override
			public void read(ByteBuffer buf, int[] vals) {
				for (int i=0; i<vals.length; i++) {
					vals[i] = buf.getInt();
				}
			}
		};
		
		public final int numBytes;
		public final int maxVal;
		
		private Encoding(int numBytes, int maxVal) {
			this.numBytes = numBytes;
			this.maxVal = maxVal;
		}
		
		public abstract void write(int[] vals, ByteBuffer buf);
		public abstract void read(ByteBuffer buf, int[] vals);
		
		public static Encoding pickBest(int maxVal) {
			for (Encoding encoding : Encoding.values()) {
				if (maxVal <= encoding.maxVal) {
					return encoding;
				}
			}
			
			// the compiler isn't smart enough to figure out that this can't happen
			throw new Error("unpossible!");
		}
	}
	
	public final RCs rcs;
	public final Encoding encoding;
	public final EntrySize entrySize;

	protected AssignmentsSerializer(RCs rcs, int numBytes) {
		this.rcs = rcs;
		
		// get the most efficient encoding, based on the biggest RC number at any position
		int maxVal = 0;
		for (int pos=0; pos<rcs.getNumPos(); pos++) {
			for (int val : rcs.get(pos)) {
				maxVal = Math.max(maxVal, val);
			}
		}
		encoding = Encoding.pickBest(maxVal);
		entrySize = EntrySize.findBigEnoughSizeFor(rcs.getNumPos()*encoding.numBytes + numBytes);
	}
	
	public EntrySize getEntrySize() {
		return entrySize;
	}
	
	protected void writeAssignments(int[] assignments, ByteBuffer buf) {

		// just in case...
		if (assignments.length != rcs.getNumPos()) {
			throw new WrongNumberOfAssignmentsException(rcs, assignments);
		}

		encoding.write(assignments, buf);
	}
	
	protected void readAssignments(ByteBuffer buf, int[] assignments) {

		// just in case...
		if (assignments.length != rcs.getNumPos()) {
			throw new WrongNumberOfAssignmentsException(rcs, assignments);
		}

		encoding.read(buf, assignments);
	}
	
	protected int[] readAssignments(ByteBuffer buf) {
		int n = rcs.getNumPos();
		int[] assignments = new int[n];
		encoding.read(buf, assignments);
		return assignments;
	}
}
