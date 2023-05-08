package osprey


enum class OS(val id: String) {

	LINUX("linux"),
	OSX("osx"),
	WINDOWS("win");

	companion object {

		operator fun get(id: String?): OS? =
			values().find { it.id == id }

		/**
		 * Lookup the operating system type using JVM property os.name
		 */
		fun get(): OS {

			// NOTE: don't rely on any internal Gradle APIs here
			// (like org.gradle.internal.os.OperatingSystem)
			// since they're not guaranteed to be stable

			val name = System.getProperty("os.name")
				?.lowercase()
				?: throw Error("no operating system defined by the JVM")

			return if (name.contains("windows")) {
				WINDOWS
			} else if (name.contains("mac os x") || name.contains("darwin") || name.contains("osx")) {
				OSX
			} else if (name.contains("linux")) {
				LINUX
			} else {
				throw Error("unrecognized operating system: $name")
			}
		}
	}
}
