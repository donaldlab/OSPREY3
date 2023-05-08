package edu.duke.cs.osprey.gui.forcefield.amber

import org.apache.commons.lang3.SystemUtils


enum class OperatingSystem {

	Linux,
	Mac,
	Windows;

	companion object {

		fun get() =
			if (SystemUtils.IS_OS_LINUX) {
				Linux
			} else if (SystemUtils.IS_OS_MAC) {
				Mac
			} else if (SystemUtils.IS_OS_WINDOWS) {
				Windows
			} else {
				throw NoSuchElementException("unrecognized operating system: ${SystemUtils.OS_NAME}")
			}
	}
}

enum class Architecture(vararg val names: String) {

	X86("x86", "i386", "i486", "i586", "i686"),
	Amd64("x86_64", "amd64");

	companion object {

		fun get(): Architecture {
			val name = System.getProperty("os.arch")
			return values()
				.find { name.lowercase() in it.names }
				?: throw NoSuchElementException("unrecognizied architecture: $name")
		}
	}
}

class Platform(
	val os: OperatingSystem,
	val arch: Architecture
) {

	companion object {
		private val platform by lazy { Platform(OperatingSystem.get(), Architecture.get()) }
		fun get() = platform
	}

	fun toPathString() = "$os/$arch".lowercase()

	override fun toString() = "$os/$arch"
}
