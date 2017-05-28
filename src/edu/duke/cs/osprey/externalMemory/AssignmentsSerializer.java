package edu.duke.cs.osprey.externalMemory;

import java.nio.ByteBuffer;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.Position;
import edu.duke.cs.tpie.EntrySize;

public abstract class AssignmentsSerializer {
	
	public static enum Encoding {
		
		// NOTE: this order is important for pickBest()
		
		Byte(java.lang.Byte.BYTES, java.lang.Byte.MAX_VALUE) {
			
			@Override
			public void write(int[] vals, ByteBuffer buf) {
				for (int i=0; i<vals.length; i++) {
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
	
	public final SimpleConfSpace space;
	public final Encoding encoding;
	public final EntrySize entrySize;
	
	protected AssignmentsSerializer(SimpleConfSpace space, int numBytes) {
		this.space = space;
		
		// get the most efficient encoding, based on the max number of RCs at any position
		int numRCs = 0;
		for (Position pos : space.positions) {
			numRCs = Math.max(numRCs, pos.resConfs.size());
		}
		encoding = Encoding.pickBest(numRCs);
		entrySize = EntrySize.findBigEnoughSizeFor(space.positions.size()*encoding.numBytes + numBytes);
	}
	
	public EntrySize getEntrySize() {
		return entrySize;
	}
	
	protected void writeAssignments(int[] assignments, ByteBuffer buf) {
		encoding.write(assignments, buf);
	}
	
	protected int[] readAssignments(ByteBuffer buf) {
		int n = space.positions.size();
		int[] assignments = new int[n];
		encoding.read(buf, assignments);
		return assignments;
	}
}
