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

import java.io.InputStream;
import java.nio.ByteBuffer;


public class ByteBufferInputStream extends InputStream {

	public final ByteBuffer buf;

	public ByteBufferInputStream(ByteBuffer buf) {
		super();
		this.buf = buf;
	}

	@Override
	public int read() {
		return buf.get() & 0xff;
	}

	@Override
	public int read(byte[] b, int off, int len) {
		int size = Math.min(len, buf.remaining());
		buf.get(b, off, size);
		if (size > 0) {
			return size;
		}
		return -1;
	}
}
