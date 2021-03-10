package edu.duke.cs.osprey.molscope.tools

import cuchaz.kludge.tools.changed
import java.util.*


/**
 * Tracks a collection of things for changes over time.
 * Things are compared by identity, not value
 * The tracker is insensitive to the order of the collection.
 */
class IdentityChangeTracker<T> {

	private val knownThings = IdentityHashMap<T,Nothing>()

	var changed: Boolean = false
		private set

	fun update(things: Collection<T>) {

		changed = things.changed(knownThings.keys)

		if (changed) {

			// save the new things
			knownThings.clear()
			things.forEach { knownThings[it] = null }
		}
	}
}
