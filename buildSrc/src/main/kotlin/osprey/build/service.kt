package osprey.build

import org.gradle.api.Project
import org.gradle.api.tasks.bundling.Compression
import org.gradle.api.tasks.bundling.Tar
import org.gradle.kotlin.dsl.creating
import org.gradle.kotlin.dsl.get
import org.gradle.kotlin.dsl.getValue

import osprey.*


object BuildService {

	// NOTE: osprey-service build scripts depend on these names, so don't change them without also updating the shell scripts
	const val name = "osprey-service"

	/**
	 * Version number of the osprey service network protocol
	 *
	 * THIS IS NOT A COSMETIC VERSION NUMBER!
	 * IT HAS A STRICTLY-ENFORCED TECHNICAL MEANING!
	 *
	 * Increment this version if the protocol changes, so clients can detect if they're compatible or not.
	 *
	 * The docker container for the service supports multiple versions simultaneously.
	 * Clients will request the version of the service they understand.
	 * So older service clients in the wild can still be supported,
	 * as long as the requested service version has been built into the docker container.
	 *
	 * For instructions on building the docker container for the Osprey service, see:
	 * docs/content/contributing/service-building.md
	 *
	 * NOTE: this line parsed by src/main/docker/service/build.sh to read the current version
	 */
	const val version = "0.3"
}


fun String.isServiceRelease(): Boolean =
	startsWith(BuildService.name) && !startsWith(BuildServiceDocker.name)
	// have to check both prefixes, since they share a common prefix themselves


fun Project.makeBuildServiceTasks() {

	val jar = tasks["jar"]

	@Suppress("UNUSED_VARIABLE")
	val serviceRelease by tasks.creating(Tar::class) {
		group = "release"
		description = "build the app server runtime for this version of the osprey service"
		dependsOn("jar")

		archiveBaseName.set(BuildService.name)
		archiveVersion.set(BuildService.version)
		destinationDirectory.set(releasesDir.toFile())
		compression = Compression.BZIP2

		val dir = buildPath / "service-${BuildService.version}"
		doFirst {
			dir.recreateFolder()

			// write the run script
			val libs = ArrayList<String>().apply {
				jar.outputs.files
					.forEach { add(it.name) }
				sourceSets.main.runtimeClasspath
					.filter { it.extension == "jar" }
					.forEach { add(it.name) }
			}
			val classpath = libs.joinToClasspath { "lib/$it" }
			writeScript(
				dir, "osprey-service",
				"""
					|cd `dirname "$0"`/..
					|java -Xmx1g -cp "$classpath" edu.duke.cs.osprey.service.MainKt $@
				""".trimMargin()
			)
		}

		into("") { // root folder
			from("README.rst")
			from("LICENSE.txt")
			from("CONTRIBUTING.rst")
		}
		into("bin") {
			from(dir)
		}
		into("lib") {
			from(jar.outputs.files)
			from(sourceSets.main.runtimeClasspath.filter { it.extension == "jar" })
		}
		into("progs") {
			from(projectDir / "progs")
		}

		// cleanup
		doLast {
			dir.deleteFolder()
		}
	}
}
