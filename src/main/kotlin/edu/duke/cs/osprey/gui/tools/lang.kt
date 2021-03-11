package edu.duke.cs.osprey.gui.tools

import kotlin.random.Random


/**
 * Returns all unique pairs of things in the list
 */
fun <T> List<T>.pairs(): List<Pair<T,T>> =
	mapIndexed { i, item ->
		subList(0, i)
			.map { item to it }
	}
	.flatten()


class UnsupportedClassException(val msg: String, val obj: Any)
	: RuntimeException("$msg: ${obj::class.simpleName}")


fun Random.nextFloatIn(min: Float, max: Float): Float =
	if (min != max) {
		nextDouble(min.toDouble(), max.toDouble()).toFloat()
	} else {
		min
	}


data class IntPair(val i1: Int, val i2: Int)


class CombineCollisionException(
	val key: String,
	val val1: String,
	val val2: String
) : RuntimeException()

/**
 * Combine two maps together.
 * Duplicate entries are collapsed into a single entry.
 * It the same key maps to multiple values, an exception is thrown.
 */
fun <K,V> combineMaps(
	map1: Map<K,V>,
	preferKeys1: Set<K>,
	ignoreKeys1: Set<K>,
	map2: Map<K,V>,
	preferKeys2: Set<K>,
	ignoreKeys2: Set<K>,
	into: MutableMap<K,V> = HashMap()
): Map<K,V> {

	// add the entries from map1
	for ((key, val1) in map1) {

		if (key in ignoreKeys1) {
			continue
		}

		into[key] = val1
	}

	// add the entries from map2, check for collisions
	for ((key, val2) in map2.entries) {

		if (key in ignoreKeys2) {
			continue
		}

		val val1 = into[key]
		into[key] = if (val1 != null && val1 != val2) {

			// break the tie with the preferred keys
			val in1 = key in preferKeys1
			val in2 = key in preferKeys2
			if (in1 && !in2) {
				val1
			} else if (!in1 && in2) {
				val2
			} else if (in1 && in2) {
				throw IllegalArgumentException("key $key is preferred by both maps, can't break tie")
			} else {
				throw CombineCollisionException(
					key.toString(),
					val1.toString(),
					val2.toString()
				)
			}
		} else {
			val2
		}
	}

	return into
}
