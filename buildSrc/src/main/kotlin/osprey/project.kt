package osprey

import org.gradle.api.Project
import java.nio.file.Path


val Project.projectPath: Path get() = projectDir.toPath().toAbsolutePath()
val Project.buildPath: Path get() = buildDir.toPath().toAbsolutePath()


/**
 * Returns true iff the command exists
 */
fun Project.commandExists(name: String): Boolean {

	// check the cache first
	commandExistence[name]?.let { return it }

	val exists = when (OS.get()) {

		OS.WINDOWS -> {
			// don't know how to do this in Windows,
			// so just assume true and hope for the best
			true
		}

		else -> {
			exec {
				commandLine("which", name)
				setIgnoreExitValue(true)
			}.exitValue == 0
		}
	}

	// update the cache
	commandExistence[name] = exists

	return exists
}

val commandExistence = HashMap<String,Boolean>()

fun Project.commandExistsOrThrow(name: String) {
	if (!commandExists(name)) {
		throw Error("command not available: $name")
	}
}


fun Project.propertyOrNull(name: String): Any? =
	try {
		property(name)
	} catch (ex: groovy.lang.MissingPropertyException) {
		null
	}

