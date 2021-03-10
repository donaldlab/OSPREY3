 package edu.duke.cs.osprey.molscope.tools

import java.util.*


/**
 * A generic bijection map, where the two related sets are named A and B
 *
 * This is, a map with an equal number of A items and B items that is 1-to-1 and onto.
 */
open class Bijection<T> {

	private val a2b = IdentityHashMap<T,T>()
	private val b2a = IdentityHashMap<T,T>()

	fun add(a: T, b: T) {
		a2b[a] = b
		b2a[b] = a
	}

	fun addAll(other: Bijection<T>) {
		a2b.putAll(other.a2b)
		b2a.putAll(other.b2a)
	}

	fun getB(a: T) = a2b[a]
	fun getA(b: T) = b2a[b]

	fun getBOrThrow(a: T) =
		getB(a) ?: throw NoSuchElementException("item $a is not in the A side")
	fun getAOrThrow(b: T) =
		getA(b) ?: throw NoSuchElementException("item $b is not in the B side")

	fun removeA(a: T): Boolean {
		val b = getB(a) ?: return false
		a2b.remove(a)
		b2a.remove(b)
		return true
	}

	fun removeB(b: T): Boolean {
		val a = getA(b) ?: return false
		b2a.remove(b)
		a2b.remove(a)
		return true
	}

	fun removeAOrThrow(a: T) {
		val b = getBOrThrow(a)
		a2b.remove(a)
		b2a.remove(b)
	}

	fun removeBOrThrow(b: T) {
		val a = getAOrThrow(b)
		b2a.remove(b)
		a2b.remove(a)
	}
}
