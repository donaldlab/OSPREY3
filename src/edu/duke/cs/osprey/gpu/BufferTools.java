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

package edu.duke.cs.osprey.gpu;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.nio.LongBuffer;

public class BufferTools {

	public static class OutOfDirectMemoryError extends VirtualMachineError {

		public OutOfDirectMemoryError(int bytes, OutOfMemoryError err) {
			super("can't allocate " + bytes + " bytes of direct memory", err);
		}
	}
	
	public static enum Type {
		
		Normal {
			
			@Override
			public ByteBuffer makeByte(int size) {
				return ByteBuffer.allocate(size);
			}
			
			@Override
			public DoubleBuffer makeDouble(int size) {
				return DoubleBuffer.allocate(size);
			}
			
			@Override
			public IntBuffer makeInt(int size) {
				return IntBuffer.allocate(size);
			}
			
			@Override
			public LongBuffer makeLong(int size) {
				return LongBuffer.allocate(size);
			}
		},
		Direct {
			
			@Override
			public ByteBuffer makeByte(int size) {
				try {
					return ByteBuffer.allocateDirect(size).order(ByteOrder.nativeOrder());
				} catch (OutOfMemoryError ex) {
					throw new OutOfDirectMemoryError(size, ex);
				}
			}
			
			@Override
			public DoubleBuffer makeDouble(int size) {
				return makeByte(size*Double.BYTES).asDoubleBuffer();
			}
			
			@Override
			public IntBuffer makeInt(int size) {
				return makeByte(size*Integer.BYTES).asIntBuffer();
			}
			
			@Override
			public LongBuffer makeLong(int size) {
				return makeByte(size*Long.BYTES).asLongBuffer();
			}
		};
		
		public abstract ByteBuffer makeByte(int size);
		public abstract DoubleBuffer makeDouble(int size);
		public abstract IntBuffer makeInt(int size);
		public abstract LongBuffer makeLong(int size);
	}
}
