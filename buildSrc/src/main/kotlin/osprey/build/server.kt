package osprey.build

import org.gradle.api.Project
import org.gradle.api.tasks.Exec
import org.gradle.api.tasks.bundling.AbstractArchiveTask
import org.gradle.api.tasks.bundling.Compression
import org.gradle.api.tasks.bundling.Tar
import org.gradle.api.tasks.bundling.Zip
import org.gradle.kotlin.dsl.creating
import org.gradle.kotlin.dsl.get
import org.gradle.kotlin.dsl.getValue

import osprey.*
import java.nio.file.Paths


object BuildServer : Build {

	override val name = "osprey-server"

	override fun getRelease(filename: String): Release? {

		// filenames look like, eg:
		//   osprey-server-linux-4.0.tbz2
		//   osprey-server-osx-4.0.tbz2
		//   osprey-server-windows-4.0.zip

		val (base, _) = Paths.get(filename).baseAndExtension()
		val parts = base.split('-')

		// get the OS
		val os = OS[parts.getOrNull(2)]
			?: run {
				System.err.println("unrecognized os for server release: $filename")
				return null
			}

		// get the version
		val version = parts.getOrNull(3)
			?.let { Version.of(it) }
			?: run {
				System.err.println("unrecognized version for server release: $filename")
				return null
			}

		return Release(this, version, os, filename)
	}
}


fun Project.makeBuildServerTasks() {

	val jar = tasks["jar"]

	val pythonWheel by tasks.creating(Exec::class) {
		group = "build"
		description = "Build python wheel"
		dependsOn("runtime", "processResources")

		inputs.files(jar.outputs.files)
		outputs.dir(pythonWheelDir)

		doFirst {

			// delete old cruft
			delete {
				delete(fileTree(pythonBuildDir) {
					include("*.whl")
				})
			}
			delete {
				delete(pythonWheelDir)
			}

			// copy python sources
			copy {
				from(pythonSrcDir) {
					includeEmptyDirs = false
					include("osprey/*.py")
				}
				into(pythonWheelDir.toFile())
			}

			val wheelOspreyDir = pythonWheelDir / "osprey"

			// copy the documentation
			copy {
				from(".") {
					include("README.md")
					include("CITING_OSPREY.txt")
					include("LICENSE.txt")
				}
				into(wheelOspreyDir.toFile())
			}

			// copy setup.py, but change the rootDir
			copy {
				from(pythonSrcDir) {
					include("setup.py")
				}
				into(pythonWheelDir.toFile())
				filter { line ->
					if (line.startsWith("rootDir = ")) {
						// make sure to escape backslashes in windows paths
						"rootDir = \"${projectPath.toString().replace("\\", "\\\\")}\""
					} else {
						line
					}
				}
			}

			val libDir = wheelOspreyDir / "lib"

			// copy osprey jar
			copy {
				from(jar)
				into(libDir.toFile())
			}

			// copy java libs
			copy {
				from(sourceSets.main.runtimeClasspath.files
					.filter { it.extension == "jar" }
					// TODO: filter down to "natives" jars only for this OS
				)
				into(libDir.toFile())
			}

			// copy the jre folder
			copy {
				from(project.runtime.jreDir)
				into((wheelOspreyDir / "jre").toFile())
			}

			// copy the progs folder if building in Linux
			if (OS.get() == OS.LINUX) {
				copy {
					from(project.projectPath.resolve("progs").toFile())
					into((wheelOspreyDir / "progs").toFile())
				}
			}
		}
		workingDir = pythonWheelDir.toFile()
		commandLine(pythonCmd, "setup.py", "bdist_wheel", "--dist-dir", pythonBuildDir.toString())
	}

	fun AbstractArchiveTask.configureServerRelease() {
		dependsOn(pythonWheel)

		val dir = buildPath / "server"

 		doFirst {
			dir.recreateFolder()

			// write installation scripts
			writeScript(
				dir,
				"install",
				"""
					|python -m pip uninstall -y osprey
					|python -m pip install --user osprey --find-link=wheelhouse --pre
				""".trimMargin()
			)
			writeScript(
				dir,
				"uninstall",
				"python -m pip uninstall -y osprey"
			)
		}

		archiveBaseName.set("osprey-server-${OS.get().id}")
		archiveVersion.set(project.version.toString())
		destinationDirectory.set(releasesDir.toFile())

		into("") { // project root
			from("README.md")
			from("LICENSE.txt")

			// install scripts, but with wildcards to accomodate
			// that they might end in .bat on Windows
			// sadly, we can't generate the filename and then include them in
			// the copy spec because Gradle won't let us change the copy spec in doFirst =(
			from(dir) {
				from("install*")
				from("uninstall*")
			}
		}
		listOf(
			"python.GMEC",
			"python.KStar",
			"python.ccs",
			"python.EWAKStar",
			"python.PAStE",
			"gpu"
		).forEach {
			into("examples/$it") {
				from("examples/$it")
			}
		}
		into("wheelhouse") {
			from(pythonBuildDir) {
				include("osprey-*.whl")
			}
		}
	}

	val serverReleasePosix by tasks.creating(Tar::class) {
		group = "release"
		compression = Compression.BZIP2
		configureServerRelease()
	}

	val serverReleaseWindows by tasks.creating(Zip::class) {
		group = "release"
		configureServerRelease()
	}

	@Suppress("UNUSED_VARIABLE")
	val serverRelease by tasks.creating {
		group = "release"
		description = "build the server release of osprey"

		// use tar.bz2, unless windows, then use zip
		when (OS.get()) {
			OS.LINUX, OS.OSX -> dependsOn(serverReleasePosix)
			OS.WINDOWS -> dependsOn(serverReleaseWindows)
		}
	}
}
