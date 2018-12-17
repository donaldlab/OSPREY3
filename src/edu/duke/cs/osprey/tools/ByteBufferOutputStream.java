package edu.duke.cs.osprey.tools;

import java.io.OutputStream;
import java.nio.ByteBuffer;


public class ByteBufferOutputStream extends OutputStream {

	public final ByteBuffer buf;

	public ByteBufferOutputStream(ByteBuffer buf) {
		super();
		this.buf = buf;
	}

	@Override
	public void write(int b) {
		buf.put((byte)(b & 0xff));
	}

	@Override
	public void write(byte[] b, int off, int len) {
		buf.put(b, off, len);
	}
}