package edu.duke.cs.osprey.molscope.tools

import java.util.*


val Boolean.toYesNo get() = if (this) "Yes" else "No"


fun Any?.identityHashCode() = System.identityHashCode(this)



private object Assertions {
	internal val Enabled: Boolean = javaClass.desiredAssertionStatus()
}

/**
 * A lazy version of assert() that only evaluates its condition if assertions are enabled.
 */
fun assert(lazyMessage: () -> Any = { "Assertion failed" }, lazyCondition: () -> Boolean) {
	if (Assertions.Enabled) {
		if (!lazyCondition()) {
			throw AssertionError(lazyMessage())
		}
	}
}


fun <T> Collection<T>.toIdentitySet(): MutableSet<T> {
	val out = identityHashSet<T>()
	out.addAll(this)
	return out
}

fun <T> identityHashSet(): MutableSet<T> =
	Collections.newSetFromMap(IdentityHashMap<T,Boolean>())

fun <K,V> identityHashMapOf(vararg pairs: Pair<K,V>) =
	IdentityHashMap<K,V>().apply {
		putAll(pairs)
	}

fun <T> identityHashSetOf(vararg values: T) =
	identityHashSet<T>().apply {
		addAll(values)
	}

fun <T,K,V> Iterable<T>.associateIdentity(transform: (T) -> Pair<K,V>): MutableMap<K,V> {
	return associateTo(IdentityHashMap(), transform)
}

fun <T,K,V> Iterable<T>.associateIdentityIndexed(transform: (Int, T) -> Pair<K,V>): MutableMap<K,V> {
	return associateIndexedTo(IdentityHashMap(), transform)
}

fun <T, K, V, M : MutableMap<in K, in V>> Iterable<T>.associateIndexedTo(destination: M, transform: (Int, T) -> Pair<K,V>): M {
	for ((index, element) in this.withIndex()) {
		destination += transform(index, element)
	}
	return destination
}

fun <K1,K2,V> Map<K1,V>.mapKeysIdentity(transform: (Map.Entry<K1,V>) -> K2): MutableMap<K2,V> {
	return mapKeysTo(IdentityHashMap<K2,V>(), transform)
}

fun <K,V1,V2> Map<K,V1>.mapValuesIdentity(transform: (Map.Entry<K,V1>) -> V2): MutableMap<K,V2> {
	return mapValuesTo(IdentityHashMap<K,V2>(), transform)
}
