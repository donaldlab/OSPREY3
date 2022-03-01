package osprey.build

import org.gradle.api.Project
import org.gradle.api.tasks.bundling.Tar
import org.gradle.kotlin.dsl.creating
import org.gradle.kotlin.dsl.getValue

import osprey.*
import java.net.URL
import java.nio.file.Paths


object BuildServiceDocker : Build {

	// NOTE: osprey-service build scripts depend on these names, so don't change them without also updating the shell scripts
	override val name = "osprey-service-docker"

	override fun getRelease(filename: String): Release? {

		// filenames look like, eg:
		//   osprey-service-docker-0.3.tar

		val (base, _) = Paths.get(filename).baseAndExtension()
		val parts = base.split('-')

		// get the version
		val version = parts.getOrNull(3)
			?.let { Version.of(it) }
			?: run {
				System.err.println("unrecognized version for service-docker release: $filename")
				return null
			}

		return Release(this, version, OS.LINUX, filename)
	}
}

fun Project.makeBuildServiceDockerTasks() {

	/**
	 * Add an entry here for each version of the service to include in the docker container
	 */
	val oldVersions = listOf(
		"1.0"
	)

	/**
	 * always include the current version, if not already there,
	 * but as the first entry
 	 */
	val versions =
		listOf(BuildService.version) + oldVersions.filter { it != BuildService.version }

	@Suppress("UNUSED_VARIABLE")
	val serviceDockerPrep by tasks.creating {
		group = "release"
		description = "Download all versions of the service releases from the donaldlab website, for the docker build script"
		doLast {

			fun filename(version: String) =
				"${BuildService.name}-$version.tbz2"

			// dowload previous versions from the release archvie
			for (version in versions) {
				val filename = filename(version)
				val path = releasesDir / filename
				if (!path.exists()) {
					println("Downloading $filename ...")
					URL(releaseArchiveUrl, filename).download(path)
				}
			}

			// unpack them all to a folder that docker can find
			val versionsDir = buildPath / "docker" / "versions"
			versionsDir.createFolderIfNeeded()

			for (version in versions) {
				val filename = filename(version)
				val path = releasesDir / filename
				val versionDir = versionsDir / "v$version"
				copy {
					from(tarTree(path.toFile()))
					into(versionDir)
				}
			}

			// write out the versions to a file so the docker build script can read it
			(buildPath / "docker" / "versions.txt").write {
				for (version in versions) {

					// write the current version on the first line
					write("$version\n")

					// write all the supported versions on the next line
					// in a json list-like format
					write("[${versions.joinToString(",")}]")
				}
			}
		}
	}

	@Suppress("UNUSED_VARIABLE")
	val serviceDockerRelease by tasks.creating(Tar::class) {
		group = "release"
		description = "build the release of the docker image for the osprey service"

		archiveBaseName.set(BuildServiceDocker.name)
		archiveVersion.set(BuildService.version)
		destinationDirectory.set(releasesDir.toFile())

		// don't bother compressing this tar
		// the VAST MAJORITY of the space is taken up by the docker image, which is already compressed
		// we won't gain much more by compressing a few small text files
		//compression = Compression.BZIP2

		val imagePath = buildPath / "docker/osprey-service-docker-${BuildService.version}.tar.bz2"
		val serviceDir = projectPath / "buildSrc/src/main/docker/service"

		val dir = buildPath / "service-docker"
		doFirst {
			dir.recreateFolder()

			// make sure the docker image has been built
			if (!imagePath.exists()) {
				throw Error("""
					|Docker image not built yet. (expected at $imagePath)
					|Gradle can't build the Docker image because Docker requires special privileges.
					|Run the build script in $serviceDir with sudo.
				""".trimMargin())
			}
		}

		into("") { // root folder
			from("README.rst")
			from("LICENSE.txt")
			from("CONTRIBUTING.rst")
			from(serviceDir / "osprey-service")
			from(serviceDir / "install.sh")
			from(serviceDir / "uninstall.sh")
			from(imagePath)
		}

		// cleanup
		doLast {
			dir.deleteFolder()
		}
	}


	/**
	 * Dear future me:
	 *
	 * You may be tempted to put a task in here to build docker images for Osprey.
	 *
	 * Don't do it!
	 *
	 * Due to the way docker works, the build steps must be run with sudo or you will get a "Permission denied" error.
	 * Don't put your user in the `docker` group to try to avoid the sudo requirement, it's a huge security risk, see:
	 * https://fosterelli.co/privilege-escalation-via-docker.html
	 *
	 * The only other option is then to run Gradle with root access.
	 * This is a Very Bad Idea, since Gradle routinely downloads code from the internet and immediately executes it.
	 * Especially when you're running it for the first time as the root user and Gradle has to entirely repopulate
	 * its caches for a new user. You wouldn't give the internet root access to your machine. Don't run Gradle
	 * as root unless you're in a VM or container.
	 *
	 * Don't subject your development machine to huge security risks for a tiny bit of convenience.
	 * Just run your docker build steps directly in a short shell script under sudo.
	 */

	@Suppress("UNUSED_VARIABLE")
	val downloadServiceReleases by tasks.creating {
		group = "release"
		description = "Download all versions of the service releases, for the docker build script"
		doLast {
			sftp {

				// what releases do we have already?
				val localReleases = releasesDir.listFiles()
					.mapNotNull { Builds.getRelease(it.fileName.toString()) }
					.filter { it.build === Builds.service }
					.toSet()

				// what releases do we need?
				val missingReleases = ls(releaseArchiveDir.toString())
					.filter { !it.attrs.isDir }
					.mapNotNull { Builds.getRelease(it.filename) }
					.filter { it.build === Builds.service && it !in localReleases }

				// download the missing releases
				if (missingReleases.isNotEmpty()) {
					for (release in missingReleases) {
						get(
							(releaseArchiveDir / release.filename).toString(),
							(releasesDir / release.filename).toString(),
							SftpProgressLogger()
						)
					}
				} else {
					println("No extra service releases to download")
				}
			}
		}
	}
}
