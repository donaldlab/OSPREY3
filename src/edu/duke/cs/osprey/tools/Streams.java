package edu.duke.cs.osprey.tools;

import java.util.stream.Stream;

public class Streams {

	public static <T> Iterable<T> toIterable(Stream<T> stream) {
		return () -> stream.iterator();
	}
}
