package edu.duke.cs.osprey.tools;

import org.tukaani.xz.LZMA2Options;
import org.tukaani.xz.XZInputStream;
import org.tukaani.xz.XZOutputStream;

import java.io.*;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;


/**
 * A utility for LZMA2 compression and decompression.
 */
public class LZMA2 {

	private static final LZMA2Options options = new LZMA2Options();
	private static final Charset charset = StandardCharsets.UTF_8;

	public static byte[] compress(byte[] bytes) {
		ByteArrayOutputStream buf = new ByteArrayOutputStream();
		try (DataOutputStream out = new DataOutputStream(new XZOutputStream(buf, options))) {

			// encode the uncompressed length as a simple int
			out.writeInt(bytes.length);

			// write the payload
			out.write(bytes);

		} catch (IOException ex) {
			throw new RuntimeException("can't compress", ex);
		}
		return buf.toByteArray();
	}

	public static byte[] compress(String text) {
		ByteArrayOutputStream buf = new ByteArrayOutputStream();
		try (XZOutputStream out = new XZOutputStream(buf, options)) {

			// encode the text as bytes
			byte[] bytes = text.getBytes(charset);

			// encode the uncompressed length as characters,
			// so if we try to view the compressed file with standard desktop tools,
			// the decompressed blob will look like a text file
			// 12 characters should be more than enough for an int in decimal encoding
			byte[] lenbuf = String.format("%-12d", bytes.length).getBytes(charset);
			if (lenbuf.length != 12) {
				throw new Error("length encoded to unexpected size: " + lenbuf.length);
			}
			out.write(lenbuf);

			// write the payload
			out.write(bytes);

		} catch (IOException ex) {
			throw new RuntimeException("can't compress", ex);
		}
		return buf.toByteArray();
	}

	public static byte[] decompressBytes(byte[] bytes) {
		try (ByteArrayInputStream bin = new ByteArrayInputStream(bytes)) {
			try (DataInputStream in = new DataInputStream(new XZInputStream(bin))) {

				// read the length, encoded as a simple int
				int length = in.readInt();

				// read the payload
				byte[] out = new byte[length];
				in.readFully(out);

				// the docs say to try to read an extra byte, so the library can check the footers
				if (in.read() != -1) {
					throw new IOException("Expected EOF, but it wasn't");
				}

				return out;
			}
		} catch (IOException ex) {
			throw new RuntimeException("can't decompress", ex);
		}
	}

	public static String decompressString(byte[] bytes) {
		try (ByteArrayInputStream bin = new ByteArrayInputStream(bytes)) {
			try (DataInputStream in = new DataInputStream(new XZInputStream(bin))) {

				// read the length, decoding the characters, then decoding the decimal
				byte[] lenbuf = new byte[12];
				in.readFully(lenbuf);
				int length = Integer.parseInt(new String(lenbuf, charset).trim());

				// read the payload
				byte[] out = new byte[length];
				in.readFully(out);

				// the docs say to try to read an extra byte, so the library can check the footers
				if (in.read() != -1) {
					throw new IOException("Expected EOF, but it wasn't");
				}

				return new String(out, charset);
			}
		} catch (IOException ex) {
			throw new RuntimeException("can't decompress", ex);
		}
	}
}
