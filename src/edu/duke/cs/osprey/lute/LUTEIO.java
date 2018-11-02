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

package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.IntEncoding;

import java.io.*;
import java.nio.charset.StandardCharsets;


public class LUTEIO {

	public static LUTEState read(File file) {
		return read(FileTools.readFileBytes(file));
	}

	public static LUTEState read(byte[] data) {
		try (DataInputStream in = new DataInputStream(new ByteArrayInputStream(data))) {

			// check the magic number
			byte[] magic = new byte[4];
			in.read(magic);
			if (!new String(magic, StandardCharsets.US_ASCII).equals("LUTE")) {
				throw new IOException("not a LUTE file");
			}

			// read the version
			byte version = in.readByte();
			switch (version) {
				case 1: return readV1(in);
				default: throw new IOException("unrecognized version: " + version);
			}

		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
	}

	public static void write(LUTEState state, File file) {
		FileTools.writeFileBytes(write(state), file);
	}

	public static byte[] write(LUTEState state) {
		ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		try (DataOutputStream out = new DataOutputStream(bytes)) {

			// write a magic number for the file
			out.write("LUTE".getBytes(StandardCharsets.US_ASCII));

			// write the version header
			out.writeByte(1);

			writeV1(state, out);

		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
		return bytes.toByteArray();
	}

	private static LUTEState readV1(DataInput in)
	throws IOException {

		// read the encodings
		IntEncoding sizeEncoding = IntEncoding.values()[in.readByte()];
		IntEncoding posEncoding = IntEncoding.values()[in.readByte()];
		IntEncoding rcEncoding = IntEncoding.values()[in.readByte()];

		// read the tuples and energies
		int n = in.readInt();
		LUTEState state = new LUTEState(n);
		for (int i=0; i<n; i++) {

			int size = sizeEncoding.read(in);
			RCTuple tuple = new RCTuple();
			for (int j=0; j<size; j++) {
				tuple.pos.add(posEncoding.read(in));
				tuple.RCs.add(rcEncoding.read(in));
			}
			state.tuples[i] = tuple;

			state.tupleEnergies[i] = in.readDouble();
		}

		// read the offset
		state.tupleEnergyOffset = in.readDouble();

		return state;
	}

	private static void writeV1(LUTEState state, DataOutput out)
	throws IOException {

		// determine the encodings for positions and RCs
		int maxSize = 0;
		int maxPos = 0;
		int maxRC = 0;
		for (RCTuple tuple : state.tuples) {
			maxSize = Math.max(maxSize, tuple.size());
			for (int i=0; i<tuple.size(); i++) {
				maxPos = Math.max(maxPos, tuple.pos.get(i));
				maxRC = Math.max(maxRC, tuple.RCs.get(i));
			}
		}
		IntEncoding sizeEncoding = IntEncoding.get(maxSize);
		IntEncoding posEncoding = IntEncoding.get(maxPos);
		IntEncoding rcEncoding = IntEncoding.get(maxRC);
		out.writeByte(sizeEncoding.ordinal());
		out.writeByte(posEncoding.ordinal());
		out.writeByte(rcEncoding.ordinal());

		// then write all the tuples and energies
		out.writeInt(state.tuples.length);
		for (int i=0; i<state.tuples.length; i++) {

			RCTuple tuple = state.tuples[i];
			sizeEncoding.write(out, tuple.size());
			for (int j=0; j<tuple.size(); j++) {
				posEncoding.write(out, tuple.pos.get(j));
				rcEncoding.write(out, tuple.RCs.get(j));
			}

			out.writeDouble(state.tupleEnergies[i]);
		}

		// write the offset
		out.writeDouble(state.tupleEnergyOffset);
	}
}
