package edu.duke.cs.osprey.gui

import org.slf4j.LoggerFactory
import java.nio.charset.Charset
import java.util.*


object OspreyGui {

	const val name = "Osprey GUI"

	val log = LoggerFactory.getLogger(OspreyGui::class.java)

	fun getResourceAsStream(path: String) = OspreyGui.javaClass.getResourceAsStream(path)

	fun getResourceAsString(path: String, charset: Charset = Charsets.UTF_8) =
		getResourceAsStream(path).use { stream -> stream.reader(charset).readText() }

	fun getResourceAsBytes(path: String) =
		getResourceAsStream(path).use { stream -> stream.readBytes() }
}
