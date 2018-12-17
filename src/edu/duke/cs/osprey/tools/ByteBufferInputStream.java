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
