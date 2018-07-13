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
