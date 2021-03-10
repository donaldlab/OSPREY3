package edu.duke.cs.osprey.molscope.render


class HoverEffects {

	inner class Writer : AutoCloseable {

		var effect: RenderEffect? = null

		init {
			writers.add(this)
		}

		override fun close() {
			writers.remove(this)
		}
	}

	private val writers = ArrayList<Writer>()

	fun writer() = Writer()

	/**
	 * Gets the effect from the most recently-created writer.
	 */
	fun get() = writers.lastOrNull()?.effect
}
