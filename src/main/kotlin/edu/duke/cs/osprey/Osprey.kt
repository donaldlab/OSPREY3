package edu.duke.cs.osprey

import org.slf4j.LoggerFactory
import java.nio.charset.Charset
import java.util.*


object Osprey {

	const val name = "Osprey"

	private val properties =
		Properties().apply {
			getResourceAsStream("build.properties")
				?.use { load(it) }
				?: throw Error("can't find build.properties")
		}

	private fun string(name: String) = properties.getProperty(name) ?: throw NoSuchElementException("no property named $name")
	private fun bool(name: String) = string(name).toBoolean()

	val dev = bool("dev")
	val version = string("version")
	val versionService = string("versionService")

	val log = LoggerFactory.getLogger(Osprey::class.java)

	fun getResourceAsStream(path: String) = Osprey.javaClass.getResourceAsStream(path)

	fun getResourceAsString(path: String, charset: Charset = Charsets.UTF_8) =
		getResourceAsStream(path).use { stream -> stream.reader(charset).readText() }

	fun getResourceAsBytes(path: String) =
		getResourceAsStream(path).use { stream -> stream.readBytes() }
}
