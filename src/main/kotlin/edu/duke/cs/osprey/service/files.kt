package edu.duke.cs.osprey.service

import java.nio.file.Files
import java.nio.file.Path


fun String.write(path: Path) {
	path.toFile().writeText(this, Charsets.UTF_8)
}

fun Path.read(): String =
	toFile().readText(Charsets.UTF_8)

fun ByteArray.write(path: Path) {
	path.toFile().writeBytes(this)
}

fun Path.readBytes(): ByteArray =
	toFile().readBytes()


val Path.exists get() = Files.exists(this)

fun Path.deleteFolder() =
	Files.walk(this)
		.sorted(Comparator.reverseOrder())
		.forEach { Files.delete(it) }


inline fun <R> tempFolder(name: String, block: (Path) -> R): R {
	val dir = Files.createTempDirectory(name)
	try {
		return block(dir)
	} finally {
		dir.deleteFolder()
	}
}
