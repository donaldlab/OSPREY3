package osprey

import org.gradle.api.Project
import org.gradle.api.tasks.Exec
import org.gradle.kotlin.dsl.creating
import org.gradle.kotlin.dsl.get
import org.gradle.kotlin.dsl.getValue


val Project.pythonSrcDir get() = projectPath / "src/main/python"
val Project.pythonBuildDir get() = buildPath / "python"
val Project.pythonWheelDir get() = pythonBuildDir / "wheel"
val Project.pythonWheelhouseDir get() = pythonSrcDir / "wheelhouse"

/** assume python3 is available on the system, but allow overriding the actual command to invoke it */
val Project.pythonCmd get() = findProperty("OSPREY_PYTHON3")?.toString() ?: "python3"


fun Project.makePythonTasks() {

	tasks["processResources"].doLast {
		// write the version to a build file so other tools (eg python scripts) can find it
		val versionFile = buildPath / "osprey-version"
		versionFile.parent.createFolderIfNeeded()
		versionFile.write {
			write(project.version.toString())
		}
	}

	@Suppress("UNUSED_VARIABLE")
	val pythonDevelop by tasks.creating(Exec::class) {
		group = "develop"
		description = "Install python package in development mode"
		workingDir = pythonSrcDir.toFile()
		dependsOn("processResources")
		commandLine(
			pythonCmd, "-m", "pip",
			"install",
			"--user", "--editable",
			".", // path to package to install, ie osprey
			"--find-links=$pythonWheelhouseDir" // add a wheelhouse dir to find any bundled packages
		)
		doLast {
			// write the java classpath somewhere our python code can find it
			pythonBuildDir.createFolderIfNeeded()
			val classpathPath = pythonBuildDir / "classpath.txt"
			classpathPath.writeLines(sourceSets.main.runtimeClasspath.files.map { it.absolutePath })
			println("Installed development Osprey package")
		}
	}

	@Suppress("UNUSED_VARIABLE")
	val pythonUndevelop by tasks.creating(Exec::class) {
		group = "develop"
		description = "Uninstall development mode python package"
		workingDir = pythonSrcDir.toFile()
		commandLine(
			pythonCmd, "-m", "pip",
			"uninstall",
			"--yes", "osprey"
		)
		doLast {
			println("Uninstalled development Osprey package")
		}
	}
}
