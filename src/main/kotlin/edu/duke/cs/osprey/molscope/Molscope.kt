package edu.duke.cs.osprey.molscope

import java.util.*


object Molscope {

	const val name = "Molscope"

	private val properties =
		Properties().apply {
			Molscope.javaClass.getResourceAsStream("build.properties")
				?.use { load(it) }
				?: throw Error("can't find build.properties")
		}

	private fun string(name: String) = properties.getProperty(name) ?: throw NoSuchElementException("no property named $name")
	private fun bool(name: String) = string(name).toBoolean()

	val dev = bool("dev")
}
