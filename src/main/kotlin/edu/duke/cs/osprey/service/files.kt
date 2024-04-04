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


/**
 * @param name the name of the temporary directory to create
 * @param block the work to be done in the newly created directory, taken as a Path parameter
 * @return whatever is returned by the function passed as 'block', eg a result structure
 */
inline fun <R> tempFolder(name: String, block: (Path) -> R): R {
	val dir = Files.createTempDirectory(name)
	try {
		return block(dir)
	} finally {
		dir.deleteFolder()
	}
}
