package edu.duke.cs.osprey.molscope.shaders

import cuchaz.kludge.tools.toByteBuffer
import java.nio.ByteBuffer


object Shaders {

	/* NOTE:
		this class should be in the "shaders" package,
		so it can reference the shader files on the classpath
		(which are also in the "shaders" package) without
		dealing with relative/absolute resource path issues
	*/

	operator fun get(name: String): ByteBuffer {
		val url = javaClass.getResource("$name.spv")
		url?.openStream()
			?.use {
				return it.readBytes().toByteBuffer()
			}
			?: throw NoSuchElementException("didn't find shader named $name")
	}
}
