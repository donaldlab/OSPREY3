package osprey.build

import org.gradle.api.Project
import org.gradle.kotlin.dsl.creating
import org.gradle.kotlin.dsl.getValue

import osprey.*
import java.nio.file.Paths


val Project.releasesDir get() = buildPath / "releases"

/**
 * Folder (in the dlab file system) where build artifacts are saved forever.
 *
 * This folder is served on the public web at:
 * https://www2.cs.duke.edu/donaldlab/software/osprey/releases/
 *
 * So files here can be downloaded by users anywhere in the world.
 */
val releaseArchiveDir = Paths.get("/usr/project/dlab/www/donaldlab/software/osprey/releases")



object Builds {
	val desktop = BuildDesktop
	val server = BuildServer
	val service = BuildService
	val serviceDocker = BuildServiceDocker
}


fun Project.makeBuildTasks() {

	// options for the Java runtime that gets bundled into the server and desktop releases, see:
	// https://badass-runtime-plugin.beryx.org/releases/latest/
	runtime {

		options.addAll(
			"--strip-debug",
			"--compress", "2",
			"--no-header-files",
			"--no-man-pages"
		)
		modules.addAll(
			// TODO: do we really need all of these?
			"java.desktop",
			"java.xml",
			"jdk.unsupported",
			"java.logging",
			"java.sql",
			"java.naming",
			"java.management",
			"jdk.httpserver",
			"jdk.zipfs", // needed to provide jar:// file system
			"jdk.incubator.foreign" // needed for foreign memory access API
		)
	}

	@Suppress("UNUSED_VARIABLE")
	val archiveReleases by tasks.creating {
		group = "release"
		description = "Upload release builds to the dlab archive, where they are downloadable by the public"
		doLast {
			sftp {

				// get the current releases
				val archivedReleases = ls(releaseArchiveDir.toString())
					.filter { !it.attrs.isDir }
					.map { it.filename }

				// diff against the local releases
				val missingReleases = releasesDir.listFiles()
					.filter { it.fileName.toString() !in archivedReleases }
					.toList()

				if (missingReleases.isNotEmpty()) {
					for (release in missingReleases) {
						put(
							release.toString(),
							(releaseArchiveDir / release.fileName).toString(),
							SftpProgressLogger()
						)
					}
				} else {
					println("No new releases to upload")
				}
			}
		}
	}
}
