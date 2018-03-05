package edu.duke.cs.osprey.tools;

import java.util.Iterator;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class Streams {

	public static <T> Stream<T> of(Iterator<T> iter) {
		return of(() -> iter);
	}

	public static <T> Stream<T> of(Iterable<T> things) {
		return StreamSupport.stream(things.spliterator(), false);
	}
}
