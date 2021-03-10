package edu.duke.cs.osprey.molscope.gui.features


class FeatureId(id: String) {

	private fun String.normalize() = toLowerCase()

	private val str: String = id.normalize()

	override fun toString() = str
	override fun hashCode() = str.hashCode()
	override fun equals(other: Any?) = other is FeatureId && this.str == other.str
}

interface HasFeatureId {
	val id: FeatureId
}
