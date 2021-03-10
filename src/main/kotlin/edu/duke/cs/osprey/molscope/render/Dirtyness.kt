package edu.duke.cs.osprey.molscope.render


/**
 * Keeps track of "dirty" status by monitoring arbitrary state
 */
class Dirtyness {

	private var state: List<Any?>? = null

	var isDirty = true
		private set

	/**
	 * returns true if the args don't match the previous args
	 * always returns true on the first invocation
	 */
	fun update(newState: List<Any?>) {
		isDirty = state != newState
		state = newState
	}

	fun update(vararg args: Any?) = update(args.toList())
}
