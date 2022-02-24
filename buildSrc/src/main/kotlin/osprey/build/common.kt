package osprey.build

import org.gradle.api.Project
import org.gradle.kotlin.dsl.creating
import org.gradle.kotlin.dsl.getValue

import osprey.*
import java.net.URL
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

/**
 * The publicly-accessible URL for the release archive
 */
val releaseArchiveUrl = URL("https://www.cs.duke.edu/donaldlab/software/osprey/releases/")


interface Build {

	val name: String

	fun isBuild(filename: String): Boolean =
		filename.startsWith(name)

	fun getRelease(filename: String): Release?
}

object Builds {

	val desktop = BuildDesktop
	val server = BuildServer
	val service = BuildService
	val serviceDocker = BuildServiceDocker

	val all = listOf(desktop, server, service, serviceDocker)

	operator fun get(filename: String): Build? =
		all.find { it.isBuild(filename) }

	fun getRelease(filename: String): Release? =
		get(filename)?.getRelease(filename)
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


data class Release(
	val build: Build,
	val version: Version,
	val os: OS,
	val filename: String
)

data class Version(
	val major: Int,
	val minor: Int,
	val patch: Int
) : Comparable<Version> {

	companion object {

		fun of(str: String): Version {

			var major = 0
			var minor = 0
			var patch = 0

			val parts = str.split(".")
			if (parts.size > 0) {
				major = parts[0].toIntOrNull()
					?: throw Error("can't parse version: $str")
			}
			if (parts.size > 1) {
				minor = parts[1].toIntOrNull()
					?: throw Error("can't parse version: $str")
			}
			if (parts.size > 2) {
				patch = parts[2].toIntOrNull()
					?: throw Error("can't parse version: $str")
			}

			return Version(major, minor, patch)
		}
	}

	override fun compareTo(other: Version): Int {

		// simple lexicographical ordering: major, minor, then patch

		this.major.compareTo(other.major)
			.takeIf { it != 0 }
			?.let { return it }

		this.minor.compareTo(other.minor)
			.takeIf { it != 0 }
			?.let { return it }

		return this.patch.compareTo(other.patch)
	}
}

fun Project.analyzeReleases(): List<Release> {

	// first, collect the releases that are currently available for public download
	val archivedReleases = sftp {
		ls(releaseArchiveDir.toString())
			.filter { !it.attrs.isDir }
			.map { it.filename }
	}

	/* DEBUG: to get local filenames for testing instead of remote filenames from SSH
	val archivedReleases = releasesDir.listFiles()
		.map { it.fileName.toString() }
		.toList()
	*/

	return archivedReleases.mapNotNull {
		Builds.getRelease(it)
	}
}
