package edu.duke.cs.osprey.molscope.render

import cuchaz.kludge.tools.AutoCloser


class OffscreenRenderer(
	val vk: VulkanDevice
) : AutoCloseable {

	private val closer = AutoCloser()
	private fun <R:AutoCloseable> R.autoClose() = also { closer.add(this@autoClose) }
	override fun close() = closer.close()

	// TODO
}