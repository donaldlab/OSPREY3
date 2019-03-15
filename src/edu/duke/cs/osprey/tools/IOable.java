package edu.duke.cs.osprey.tools;

import java.io.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;


public interface IOable {

	void writeTo(DataOutput out) throws IOException;
	void readFrom(DataInput in) throws IOException;

	default void writeTo(File file) {
		try (FileOutputStream out = new FileOutputStream(file)) {
			writeTo(new DataOutputStream(out));
		} catch (IOException ex) {
			throw new RuntimeException("can't write to file: " + file.getAbsolutePath(), ex);
		}
	}

	default void readFrom(File file) {
		try (FileInputStream in = new FileInputStream(file)) {
			readFrom(new DataInputStream(in));
		} catch (IOException ex) {
			throw new RuntimeException("can't read from file: " + file.getAbsolutePath(), ex);
		}
	}

	default void writeToGZIP(File file) {
		try (FileOutputStream out = new FileOutputStream(file)) {
			GZIPOutputStream gout = new GZIPOutputStream(out);
			writeTo(new DataOutputStream(gout));
			gout.finish();
		} catch (IOException ex) {
			throw new RuntimeException("can't write to file: " + file.getAbsolutePath(), ex);
		}
	}

	default void readFromGZIP(File file) {
		try (FileInputStream in = new FileInputStream(file)) {
			GZIPInputStream gin = new GZIPInputStream(in);
			readFrom(new DataInputStream(gin));
		} catch (IOException ex) {
			throw new RuntimeException("can't read from file: " + file.getAbsolutePath(), ex);
		}
	}

	static IOable of(Object obj) {
		if (obj instanceof IOable) {
			return (IOable)obj;
		} else {
			throw new IllegalArgumentException("object " + obj.getClass().getSimpleName() + " is not IOable");
		}
	}
}
