package osprey.build

import org.gradle.api.Project
import org.gradle.api.tasks.bundling.Tar
import org.gradle.kotlin.dsl.creating
import org.gradle.kotlin.dsl.getValue

import osprey.*


object BuildServiceDocker {
	// NOTE: osprey-service build scripts depend on these names, so don't change them without also updating the shell scripts
	const val name = "osprey-service-docker"
}

fun Project.makeBuildServiceDockerTasks() {

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
				val localReleaseNames = releasesDir.listFiles()
					.map { it.fileName.toString() }
					.filter { it.isServiceRelease() }
					.toSet()

				// what releases do we need?
				val missingReleases = ls(releaseArchiveDir.toString())
					.filter { !it.attrs.isDir && it.filename.isServiceRelease() && it.filename !in localReleaseNames }

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
