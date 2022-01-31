package osprey

import java.io.*
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.PosixFilePermission
import java.util.stream.Stream
import kotlin.streams.asSequence


// some conveniences for files and paths

fun Path.exists(): Boolean =
	Files.exists(this)

fun Path.deleteFolder() =
	toFile().deleteRecursively()

fun Path.createFolderIfNeeded(): Path {
	Files.createDirectories(this)
	return this
}

fun Path.recreateFolder() {
	deleteFolder()
	createFolderIfNeeded()
}

fun Path.copyTo(to: Path) =
	Files.copy(this, to)

fun Path.copyFolderTo(to: Path) =
	toFile().copyRecursively(to.toFile(), overwrite = true)

fun Path.listFiles(): Sequence<Path> =
	Files.list(this).asSequence()
		.filter { Files.isRegularFile(it) }

fun Path.read(): String =
	Files.readString(this)

fun Path.write(block: Writer.(OutputStream) -> Unit) {
	toFile().outputStream().use { out ->
		out.writer().use { writer ->
			writer.block(out)
		}
	}
}

fun Path.writeLines(lines: Iterable<String>) {
	write {
		for (line in lines) {
			write(line)
			write("\n")
		}
	}
}

fun Path.rename(dst: String) {
	if ('/' in dst) {
		throw IllegalArgumentException("invalid rename, can't have /")
	}
	Files.move(this, parent.resolve(dst))
}

fun Path.extension(): String? =
	fileName.toString().split('.')
		.takeIf { it.size > 1 }
		?.last()

fun Path.deleteFile() =
	Files.delete(this)


operator fun Path.div(suffix: String): Path =
	resolve(suffix)

operator fun Path.div(suffix: Path): Path =
	resolve(suffix)

operator fun File.div(suffix: String): File =
	resolve(suffix)


fun <T> Iterable<T>.joinToClasspath(transform: (T) -> CharSequence = { it.toString() }): String =
	joinToString(File.pathSeparator, transform = transform)


fun writeScript(dir: Path, filename: String, cmd: String): Path =
	when (OS.get()) {
		OS.OSX,
		OS.LINUX -> {

			val file = dir / "$filename.sh"
			file.write {
				write("""
					|#! /bin/sh
					|$cmd
				""".trimMargin())
			}

			// set the shell script executable
			Files.setPosixFilePermissions(file, Files.getPosixFilePermissions(file).apply {
				add(PosixFilePermission.OWNER_EXECUTE)
			})

			file
		}
		OS.WINDOWS -> {

			val file = dir / "$filename.bat"
			file.write {
				write("""
					|@echo off
					|$cmd
				""".trimMargin())
			}

			file
		}
	}


fun <T> Path.walk(block: (Stream<Path>) -> T): T =
	Files.walk(this).use(block)
