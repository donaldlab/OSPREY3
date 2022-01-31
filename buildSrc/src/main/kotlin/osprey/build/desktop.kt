package osprey.build

import org.gradle.api.Project
import org.gradle.kotlin.dsl.creating
import org.gradle.kotlin.dsl.get
import org.gradle.kotlin.dsl.getValue

import osprey.*
import java.nio.file.Path
import java.nio.file.Paths


object BuildDesktop {
	const val name = "osprey-desktop"
}

fun Project.makeBuildDesktopTasks() {

	fun checkJpackage(path: Path) {
		if (!path.exists()) {
			throw NoSuchElementException("""
				|jpackage was not found at the path given: $path
				|Make sure JDK 14 (or higher) is installed and that systemProp.jpackage.home path points to its home directory.
			""".trimMargin())
		}
	}

	runtime.jpackage {

		imageName = BuildDesktop.name
		imageOutputDir = (buildPath / "jpackage").toFile()
		installerName = imageName
		mainClass = "$group.osprey.gui.MainKt"
		jpackageHome = jPackageHomeOrDefault

		val binDir = Paths.get(jpackageHome) / "bin"

		when (OS.get()) {

			OS.LINUX -> {

				checkJpackage(binDir / "jpackage")

				installerType = "deb"
				// TOOO: support `rpm` type (Fedora/RHEL/CentOS/Rocky) via command-line flag?
				// NOTE: arch linux not supported by jpackage, see:
				// https://bugs.openjdk.java.net/browse/JDK-8265498
			}

			OS.OSX -> {

				checkJpackage(binDir / "jpackage")

				installerType = "dmg"
				jvmArgs = listOf("-XstartOnFirstThread")
			}

			OS.WINDOWS -> {

				checkJpackage(binDir / "jpackage.exe")

				installerType = "msi"
				installerOptions = listOf("--win-per-user-install", "--win-dir-chooser", "--win-menu", "--win-shortcut")
				// useful for debugging launcher issues
				//imageOptions = listOf("--win-console")
			}
		}
	}


	val compileShaders by tasks.creating {
		group = "build"
		doLast {

			val outDir = sourceSets.main.resources.srcDirs.first() / "${Jvm.packagePath}/molscope/shaders"

			val inDir = file("src/main/glsl")
			fileTree(inDir)
				.matching {
					include("**/*.vert")
					include("**/*.geom")
					include("**/*.frag")
					include("**/*.comp")
				}
				.forEach { inFile ->
					val inFileRel = inFile.relativeTo(inDir)
					val outFile = inFileRel.resolveSibling(inFileRel.name + ".spv")
					exec {
						this.workingDir = outDir
						commandLine(
							"glslc",
							"--target-env=vulkan1.0",
							// some compilers don't have this flag, but it's redundant with vulkan1.0 anyway
							//"--target-spv=spv1.0",
							"-Werror",
							"-x", "glsl",
							"-o", outFile.path,
							inFile.absolutePath
						)
					}
				}
		}
	}

	tasks["build"].dependsOn(compileShaders)

	// add documentation to the gui distribution
	val desktopRelease by tasks.creating {
		group = "release"
		description = "build the desktop release of osprey"
		dependsOn(tasks.jpackage)
	}

	tasks.jpackageImage {
		doLast {
			val jp = project.runtime.jpackageData.get()
			val imageDir = when (OS.get()) {
				OS.OSX -> jp.imageOutputDir / (jp.imageName + ".app") / "Contents"
				else -> jp.imageOutputDir / jp.imageName
			}
			copy {
				from(
					projectDir / "readme.md",
					projectDir / "LICENSE.txt"
				)
				into(imageDir)
			}
		}
	}
}
