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

import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.Objects;
import java.util.function.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class Streams {

	public static <T> Stream<T> of(Iterator<T> iter) {
		return of(() -> iter);
	}

	public static <T> Stream<T> of(Iterable<T> things) {
		return StreamSupport.stream(things.spliterator(), false);
	}

	public static <T> Stream<T> of(T[] things) {
		return Arrays.stream(things);
	}

	public static <T> String joinToString(Iterable<T> things, String delimiter, Function<T,String> map) {
		if (things == null) {
			return "";
		}
		return String.join(delimiter, of(things)
			.map(map)
			.collect(Collectors.toList())
		);
	}

	public static <T> String joinToString(Iterable<T> things, String delimiter) {
		return joinToString(things, delimiter, thing -> Objects.toString(thing));
	}

	public static <T> String joinToString(T[] things, String delimiter, Function<T,String> map) {
		if (things == null) {
			return "";
		}
		return String.join(delimiter, of(things)
			.map(map)
			.collect(Collectors.toList())
		);
	}

	public static <T> String joinToString(T[] things, String delimiter) {
		return joinToString(things, delimiter, thing -> Objects.toString(thing));
	}

	public static <K,V> String joinToString(Map<K,V> things, String delimiter, BiFunction<K,V,String> map) {
		if (things == null) {
			return "";
		}
		return String.join(delimiter, of(things.entrySet())
			.map(entry -> map.apply(entry.getKey(), entry.getValue()))
			.collect(Collectors.toList())
		);
	}

	public static <K,V> String joinToString(Map<K,V> things, String delimiter) {
		return joinToString(things, delimiter, (key, val) -> Objects.toString(key) + "=" + Objects.toString(val));
	}

	public static String joinToString(int[] things, String delimiter, IntFunction<String> map) {
		if (things == null) {
			return "";
		}
		return String.join(delimiter, Arrays.stream(things)
			.mapToObj(map)
			.collect(Collectors.toList())
		);
	}

	public static String joinToString(long[] things, String delimiter, LongFunction<String> map) {
		if (things == null) {
			return "";
		}
		return String.join(delimiter, Arrays.stream(things)
			.mapToObj(map)
			.collect(Collectors.toList())
		);
	}

	public static String joinToString(double[] things, String delimiter, DoubleFunction<String> map) {
		if (things == null) {
			return "";
		}
		return String.join(delimiter, Arrays.stream(things)
			.mapToObj(map)
			.collect(Collectors.toList())
		);
	}
}
