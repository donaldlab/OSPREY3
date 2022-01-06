/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

import java.io.BufferedWriter
import java.io.OutputStreamWriter
import java.io.ByteArrayOutputStream
import java.nio.charset.StandardCharsets
import java.nio.file.attribute.PosixFilePermission
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import kotlin.streams.asSequence
import org.gradle.internal.os.OperatingSystem
import org.jetbrains.kotlin.gradle.tasks.KotlinCompile
import com.jcraft.jsch.JSch
import com.jcraft.jsch.ChannelSftp
import com.jcraft.jsch.Logger as JschLogger
import com.jcraft.jsch.SftpProgressMonitor


plugins {
	`java-library`
	kotlin("jvm") version "1.6.0"
	kotlin("plugin.serialization") version "1.6.0"
	application
	idea
	id("org.openjfx.javafxplugin") version("0.0.7")
	id("org.beryx.runtime") version "1.12.5"
}

buildscript {
	repositories {
		mavenCentral()
	}
	dependencies {
		// SSH client, BSD license: https://github.com/mwiede/jsch
		classpath("com.github.mwiede:jsch:0.1.66")
	}
}

javafx {
	version = "11"
	modules("javafx.controls")
}


val projectDir = project.projectDir.toPath().toAbsolutePath()
val buildDir = project.buildDir.toPath().toAbsolutePath()
val pythonSrcDir = projectDir.resolve("src/main/python")
val pythonBuildDir = buildDir.resolve("python")
val pythonWheelDir = pythonBuildDir.resolve("wheel")
val pythonWheelhouseDir = pythonSrcDir.resolve("wheelhouse")
val docSrcDir = pythonSrcDir.resolve("doc")
val docBuildDir = pythonBuildDir.resolve("doc")
val versionFile = buildDir.resolve("osprey-version")
val docDir = projectDir.resolve("doc")
val docMainDir = docDir.resolve("content/documentation/main")
val releasesDir = buildDir.resolve("releases")

// NOTE: osprey-service build scripts depend on these names, so don't change them without also updating the shell scripts
val releaseNameService = "osprey-service"
val releaseNameServiceDocker = "osprey-service-docker"

fun String.isServiceRelease(): Boolean =
	startsWith(releaseNameService) && !startsWith(releaseNameServiceDocker)
	// have to check both prefixes, since they share a common prefix themselves

/**
 * Folder (in the dlab file system) where build artifacts are saved forever.
 *
 * This folder is served on the public web at:
 * https://www2.cs.duke.edu/donaldlab/software/osprey/releases/
 *
 * So files here can be downloaded by users anywhere in the world.
 */
val releaseArchiveDir = Paths.get("/usr/project/dlab/www/donaldlab/software/osprey/releases")

group = "edu.duke.cs"

/**
 * Version number for Osprey itself
 *
 * This version number is largely cosmetic.
 * But it does help users provide some information to developers when reporting issues.
 */
version = "3.2"

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
val versionService = "0.3"

val packagePath = "edu/duke/cs/osprey"

// add the module dependencies directly to the javac args
// I don't think gradle has a good way to handle this yet?
val moduleArgs = listOf(
		"--add-modules=jdk.incubator.foreign"
)
// TODO: add the module args for distributions too?


repositories {
	mavenCentral()
}

val javaLangVersion = 17
java {
	toolchain {
		languageVersion.set(JavaLanguageVersion.of(javaLangVersion))
	}
}

idea {
	module {
		// use the same output folders as gradle, so the pythonDevelop task works correctly
		outputDir = sourceSets["main"].output.classesDirs.files.first()
		testOutputDir = sourceSets["test"].output.classesDirs.files.first()
		inheritOutputDirs = false
	}
}

application {
	mainClassName = "edu.duke.cs.osprey.design.Main"
}

dependencies {

	// kotlin runtime
	implementation(kotlin("stdlib-jdk8"))
	implementation(kotlin("reflect"))

	// test dependencies
	testImplementation("org.hamcrest:hamcrest-all:1.3")
	testImplementation("junit:junit:4.12")
	testImplementation("io.kotlintest:kotlintest-runner-junit5:3.4.0")

	// handle logging
	implementation("ch.qos.logback:logback-classic:1.2.3")
	implementation("org.slf4j:jul-to-slf4j:1.7.30")

	// internal osprey libs
	implementation("colt:colt:1.2.0")
	implementation("org.apache.commons:commons-math3:3.6.1")
	implementation("org.apache.commons:commons-collections4:4.1")
	implementation("com.joptimizer:joptimizer:3.5.1")
	implementation("org.ojalgo:ojalgo:41.0.0")
	implementation("org.jogamp.gluegen:gluegen-rt:2.3.2")
	implementation("org.jogamp.jocl:jocl:2.3.2")
	implementation("org.mapdb:mapdb:3.0.5")
	implementation("org.apache.xmlgraphics:batik-svggen:1.9.1")
	implementation("org.apache.xmlgraphics:batik-svg-dom:1.9.1")
	implementation("com.github.haifengl:smile-core:1.5.1")
	implementation("com.github.haifengl:smile-netlib:1.5.1")
	implementation("de.lmu.ifi.dbs.elki:elki:0.7.1")
	implementation("ch.obermuhlner:big-math:2.3.0")
	implementation("org.joml:joml:1.9.19")
	implementation("org.tukaani:xz:1.8")
	implementation("com.hazelcast:hazelcast:4.0")
	implementation("net.java.dev.jna:jna:5.5.0")
	implementation("com.google.guava:guava:29.0-jre")
	implementation("org.apache.commons:commons-lang3:3.4")
	implementation("commons-io:commons-io:2.5")
	implementation("org.tomlj:tomlj:1.0.0")
	implementation(fileTree("lib"))

	val autoValueVersion = 1.7
	implementation("ch.obermuhlner:big-math:2.3.0")
	implementation("com.fasterxml.jackson.core:jackson-databind:2.9.9")
	implementation("com.fasterxml.jackson.dataformat:jackson-dataformat-yaml:2.9.9")
	implementation("com.beust:jcommander:1.72")
	implementation("one.util:streamex:0.7.3")
	implementation(platform("software.amazon.awssdk:bom:2.15.48"))
	implementation("software.amazon.awssdk:s3")
	implementation("org.postgresql:postgresql:42.2.16")
	implementation("org.sql2o:sql2o:1.6.0")
	implementation("com.google.auto.value:auto-value-annotations:$autoValueVersion")
	annotationProcessor("com.google.auto.value:auto-value:$autoValueVersion")
	testImplementation("org.junit.jupiter:junit-jupiter:5.4.2")
	testImplementation("org.assertj:assertj-core:3.18.1")


	val ktorVersion = "1.5.4"

	// used by the gui
	implementation("com.cuchazinteractive:kludge:0.2")
	implementation("io.ktor:ktor-client-apache:$ktorVersion")
	implementation("io.ktor:ktor-client-serialization-jvm:$ktorVersion")

	// used by the service
	implementation("io.ktor:ktor-server-netty:$ktorVersion")
	implementation("io.ktor:ktor-serialization:$ktorVersion")

	// for JCuda, gradle tries (and fails) download the natives jars automatically,
	// so turn off transitive dependencies. we'll deal with natives manually
	implementation("org.jcuda:jcuda:0.8.0") {
		isTransitive = false
	}
	implementation("one.util:streamex:0.7.3")

	// native libs for GPU stuff
	listOf("natives-linux-amd64", "natives-macosx-universal", "natives-windows-amd64").forEach {
		runtimeOnly("org.jogamp.gluegen:gluegen-rt:2.3.2:$it")
		runtimeOnly("org.jogamp.jocl:jocl:2.3.2:$it")
	}
	listOf("linux-x86_64", "apple-x86_64", "windows-x86_64").forEach {
		runtimeOnly("org.jcuda:jcuda-natives:0.8.0:$it")
	}

	// used by the build system
	testImplementation("org.json:json:20210307")
	val dokkaVersion = "1.5.31"
	testImplementation("org.jetbrains.dokka:dokka-cli:$dokkaVersion")
	testImplementation("org.jetbrains.dokka:dokka-base:$dokkaVersion")
}


/* NOTE: the IDE thinks `jvmArgs` and `args` are not nullable and shows warnings
	(and the Kotlin language rules agree with that, as far as I can tell),
	but for some reason, the Kotlin compiler thinks they are nullable
	so we need the not null assertions (ie !!). Maybe a compiler bug?
*/
@Suppress("UNNECESSARY_NOT_NULL_ASSERTION")
tasks.withType<JavaExec> {

	// add the module args, if not already there
	for (arg in moduleArgs) {
		if (arg !in jvmArgs!!) {
			jvmArgs = jvmArgs!! + listOf(arg)
			// NOTE: for some reason, .add() and += don't work here
		}
	}

	doFirst {

		// make Gradle dump the java command to the console when running
		// since IntelliJ doesn't do it when using Gradle as the runner
		val jvmArgsStr = jvmArgs!!.joinToString(" ")
		val classpathStr = classpath.joinToClasspath()
		val argsStr = args!!.joinToString(" ")
		println("""
			|Java:
			|	$workingDir
			|	$executable
			|	$jvmArgsStr -cp "$classpathStr" $main $argsStr
		""".trimMargin())
	}
}


tasks.withType<JavaCompile> {
	options.compilerArgs.addAll(moduleArgs)
}

tasks.withType<KotlinCompile> {

	kotlinOptions {

		jvmTarget = javaLangVersion.toString()

		// enable experimental features so we can use the fancy ktor stuff
		freeCompilerArgs += "-Xuse-experimental=kotlin.Experimental"
		freeCompilerArgs += "-XXLanguage:+InlineClasses"
	}
}

tasks.withType<Test> {
	// the default 512m is too little memory to run test designs
	maxHeapSize = "2g"
	useJUnit()
    failFast = true

	// add the module args, if not already there
	for (arg in moduleArgs) {
		if (arg !in jvmArgs!!) {
			jvmArgs = jvmArgs!! + listOf(arg)
			// NOTE: for some reason, .add() and += don't work here
		}
	}

	testLogging {
		setExceptionFormat("full")
        events("passed", "skipped", "failed")
	}
}


// get the OS we're building on
val os = OperatingSystem.current()
val osname: String by lazy {
	when (os) {
		OperatingSystem.LINUX -> "linux"
		OperatingSystem.WINDOWS -> "win"
		OperatingSystem.MAC_OS -> "osx"
		else -> throw Error("unrecognized operating system: $os")
	}
}


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
		"jdk.incubator.foreign"
	)

	fun checkJpackage(path: Path) {
		if (!Files.exists(path)) {
			throw NoSuchElementException("""
				|jpackage was not found at the path given: $path
				|Make sure JDK 14 (or higher) is installed and that systemProp.jpackage.home path points to its home directory.
			""".trimMargin())
		}
	}

	jpackage {

		imageName = "Osprey"
		imageOutputDir = buildDir.resolve("jpackage").toFile()
		installerName = imageName
		mainClass = "$group.osprey.gui.MainKt"
		jpackageHome = jPackageHomeOrDefault

		when (os) {

			OperatingSystem.LINUX -> {

				checkJpackage(Paths.get(jpackageHome).resolve("bin").resolve("jpackage"))

				installerType = "deb"
			}

			OperatingSystem.MAC_OS -> {

				checkJpackage(Paths.get(jpackageHome).resolve("bin").resolve("jpackage"))

				installerType = "dmg"
				jvmArgs = listOf("-XstartOnFirstThread")
			}

			OperatingSystem.WINDOWS -> {

				checkJpackage(Paths.get(jpackageHome).resolve("bin").resolve("jpackage.exe"))

				installerType = "msi"
				installerOptions = listOf("--win-per-user-install", "--win-dir-chooser", "--win-menu", "--win-shortcut")
				// useful for debugging launcher issues
				//imageOptions = listOf("--win-console")
			}

			else -> throw Error("unrecognized operating system: $os")
		}
	}
}


fun isCommand(cmd: String) =
	exec {
		val output = ByteArrayOutputStream()
		isIgnoreExitValue = true
		standardOutput = output
		errorOutput = output

		when (os) {
			OperatingSystem.MAC_OS,
			OperatingSystem.LINUX -> commandLine("which", cmd)
			OperatingSystem.WINDOWS -> commandLine("powershell", "get-command", cmd)
			else -> throw Error("unrecognized operating system: $os")
		}
	}.exitValue == 0

class Python(val cmd: String) {

	// find out the version
	val version = if (isCommand(cmd)) {
		ByteArrayOutputStream().use { stdout ->
			exec {
				commandLine(cmd, "--version")
				standardOutput = stdout
				errorOutput = stdout
				isIgnoreExitValue = true
			}
			// should return something like "Python 2.7.17"
			// Except: Microsoft (in their infinite wisdom) has decided to install a fake python by default
			// This fake python writes some garbage message to stdout, so make sure we ignore it by checking the version
			stdout.toString()
				.split(" ").getOrNull(1)
				?.split(".")?.getOrNull(0)
				?.toIntOrNull()
		}
	} else {
		null
	}

	override fun toString() = "Python $version"
}

// find out what pythons are available
val pythons by lazy {
	val python3 = findProperty("OSPREY_PYTHON3")?.toString() ?: "python3"
	val python2 = findProperty("OSPREY_PYTHON2")?.toString() ?: "python2"

	listOf(
		Python(python3),
		Python(python2),
		Python("python")
	)
	.filter { it.version != null }
}

// try to find python 2 and 3
val python2 by lazy {
	pythons.find { it.version == 2 }
}
val python3 by lazy {
	pythons.find { it.version == 3 }
}

// get the default python, prefer v3 over v2
val defaultPython by lazy {
	python3
		?: python2
		?: throw Error("No python detected")
}

// only log this message when passed --info (or --debug, etc.) flag
logger.info("""
	|  Python versions found:
	|         Pythons:  ${pythons.map { it.cmd }}
	|        Python 2:  ${if (python2 != null) "✓" else "✗"}
	|        Python 3:  ${if (python3 != null) "✓" else "✗"}
	|  default Python:  $defaultPython = ${defaultPython.cmd}
""".trimMargin())


distributions {

	get("main").apply {
		distributionBaseName.set("osprey-cli")
		contents {
			into("") { // project root
				from("README.rst")
				from("LICENSE.txt")
				from("CONTRIBUTING.rst")
			}
			listOf("1CC8", "1FSV", "2O9S_L2", "2RL0.kstar", "3K75.3LQC", "4HEM", "4NPD").forEach {
				into("examples/$it") {
					from("examples/$it")
				}
			}
		}
	}

	// make distributions for both python versions
	for (version in listOf(2, 3)) {

		create("python$version").apply {
			distributionBaseName.set("osprey-$osname-python$version")
			contents {
				into("") { // project root
					from("README.rst")
					from("LICENSE.txt")
					from("CONTRIBUTING.rst")
					from(pythonBuildDir) {
						include("install.sh")
						include("install.bat")
						include("uninstall.sh")
						include("uninstall.bat")
					}
				}
				into("doc") {
					from(docBuildDir) {
						exclude(".doctrees")
					}
				}
				listOf("python.GMEC", "python.KStar", "gpu").forEach {
					into("examples/$it") {
						from("examples/$it")
					}
				}
				into("wheelhouse") {
					if (version == 2) {
						from(pythonWheelhouseDir) {
							when (os) {
								OperatingSystem.MAC_OS -> include("JPype_py2-*-macosx_*.whl")
								OperatingSystem.LINUX -> include("JPype_py2-*-linux_*.whl")
								OperatingSystem.WINDOWS -> include("JPype_py2-*-win_*.whl")
							}
						}
					}
					from(pythonBuildDir) {
						include("osprey-*.whl")
					}
				}
			}
		}
	}

	// for running tests on servers
	create("test").apply {
		distributionBaseName.set("osprey-test")
		contents {

			into("") { // project root
				from("LICENSE.txt")
			}

			// include libs, classes, and resources
			val files = sourceSets["test"].runtimeClasspath
			into("lib") {
				from(files
					.filter { it.extension == "jar" }
				)
			}
			into("classes") {
				from(files
					.filter { it.isDirectory }
					// exclude resources dirs, they're apparently already in the classes dirs
					// except this time they're not?
					//.filter { !it.endsWith("resources/main") }
					//.filter { !it.endsWith("resources/test") }
				)
			}

			// add the run script
			into("bin") {
				from(buildDir) {
					include("run.*")
				}
			}
		}
	}
}


// assume we're doing a dev build if the top task is "classes"
if (gradle.startParameter.taskNames.any { it.endsWith(":classes") }) {
	System.setProperty("isDev", true.toString())
}
val isDev = object {
	override fun toString() = System.getProperty("isDev") ?: false.toString()
}


tasks {

	// turn off tar files for all distributions
	filterIsInstance<Tar>()
		.forEach { it.enabled = false }

	val testClasses = "testClasses" {}

	processResources {

		// always update the build properties, versions, etc
		outputs.upToDateWhen { false }

		var version = version.toString()

		// append the CI build ID, if available
		version += if (rootProject.hasProperty("AZURE_BUILD_ID")) {
			val versionId = rootProject.property("AZURE_BUILD_ID")
			".$versionId"
		} else {
			"-dev"
		}

		// write the build properties so the app can find them
		from(sourceSets["main"].resources.srcDirs) {
			include("$packagePath/build.properties")
			expand(
				"dev" to isDev,
				"version" to version,
				"versionService" to versionService
			)
		}

		// I have no idea why Gradle thinks the build properties are duplicated
		// (there's only one instance of build.properties anywhere in the project),
		// but make sure we include the expanded version of the file rather than ignoring it
		duplicatesStrategy = DuplicatesStrategy.INCLUDE

		doLast {
			// write the version to a build file so other tools (eg python scripts) can find it
			versionFile.parent.toFile().mkdirs()
			versionFile.toFile().writeText(version)
		}
	}

	val jar = "jar" {}.get()

	val compileCuda_residueForcefield by creating(Exec::class) {
		nvcc(this, "residueForcefield")
	}

	val compileCuda_residueCcd by creating(Exec::class) {
		nvcc(this, "residueCcd", maxRegisters=64)
	}

	val compileCuda by creating {
		description = "Compile cuda kernels"
		dependsOn(
			compileCuda_residueForcefield,
			compileCuda_residueCcd
		)
	}

	val cleanDoc by creating(Delete::class) {
		group = "documentation"
		description = "Cleans python documentation"
		delete(docBuildDir)
	}
	val makeDoc by creating(Exec::class) {
		group = "documentation"
		description = "Build python documentation"
		workingDir = docSrcDir.toFile()
		dependsOn(processResources)
		commandLine("sphinx-build", "-b", "html", ".", "$docBuildDir")
	}
	val remakeDoc by creating {
		group = "documentation"
		description = "runs cleanDoc, then makeDoc to refresh all changes to documentation"
		dependsOn(cleanDoc, makeDoc)
	}

	val pythonDevelop by creating(Exec::class) {
		group = "develop"
		description = "Install python package in development mode"
		workingDir = pythonSrcDir.toFile()
		dependsOn(processResources)
		commandLine(
			defaultPython.cmd, "-m", "pip",
			"install",
			"--user", "--editable",
			".", // path to package to install, ie osprey
			"--find-links=$pythonWheelhouseDir" // add a wheelhouse dir to find our Jpype-py2
		)
		doLast {
			Files.createDirectories(pythonBuildDir)
			val classpathPath = pythonBuildDir.resolve("classpath.txt")
			Files.write(classpathPath, sourceSets["main"].runtimeClasspath.files.map { it.toString() })
			println("Installed development Osprey package for $defaultPython")
		}
	}

	val pythonUndevelop by creating(Exec::class) {
		group = "develop"
		description = "Uninstall development mode python package"
		workingDir = pythonSrcDir.toFile()
		commandLine(
			defaultPython.cmd, "-m", "pip",
			"uninstall",
			"--yes", "osprey"
		)
		doLast {
			println("Uninstalled development Osprey package for $defaultPython")
		}
	}

	fun Exec.pythonWheel(python: Python) {
		group = "build"
		description = "Build $python wheel"
		dependsOn("runtime", processResources)

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

			val wheelOspreyDir = pythonWheelDir.resolve("osprey")

			// copy the documentation
			copy {
				from(".") {
					include("*.rst")
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
						"rootDir = \"${projectDir.toString().replace("\\", "\\\\")}\""
					} else {
						line
					}
				}
			}

			val libDir = wheelOspreyDir.resolve("lib")

			// copy osprey jar
			copy {
				from(jar)
				into(libDir.toFile())
			}

			// copy java libs
			copy {
				from(sourceSets["main"].runtimeClasspath.files
					.filter { it.extension == "jar" }
					// TODO: filter down to "natives" jars only for this OS
				)
				into(libDir.toFile())
			}

			// copy the jre folder
			copy {
				from(runtime.get().jreDir)
				into(wheelOspreyDir.resolve("jre").toFile())
			}
		}
		workingDir = pythonWheelDir.toFile()
		commandLine(python.cmd, "setup.py", "bdist_wheel", "--dist-dir", pythonBuildDir.toString())
	}

	val python2Wheel by creating(Exec::class) {
		python2
			?.let { pythonWheel(it) }
			?: doLast {
				throw Error("no python 2 detected")
			}
	}
	val python3Wheel by creating(Exec::class) {
		python3
			?.let { pythonWheel(it) }
			?: doLast {
				throw Error("no python 3 detected")
			}
	}

	val python2InstallScripts by creating {
		group = "build"
		description = "Make install scripts for python 2 distribution"
		doLast {
			writeScript(
				pythonBuildDir, "install",
				"""
					|python -m pip uninstall -y osprey JPype-py2
					|python -m pip install --user 'numpy>=1.6,<1.16'
					|python -m pip install --user osprey --no-index --find-link=wheelhouse --pre
				""".trimMargin()
			)
		}
	}

	val python3InstallScripts by creating {
		group = "build"
		description = "Make install scripts for python 3 distribution"
		doLast {
			writeScript(
				pythonBuildDir, "install",
				"""
					|python -m pip uninstall -y osprey
					|python -m pip install --user osprey --find-link=wheelhouse --pre
				""".trimMargin()
			)
		}
	}

	val python2UninstallScripts by creating {
		group = "build"
		description = "Make uninstall scripts for python 2 distribution"
		doLast {
			writeScript(
				pythonBuildDir, "uninstall",
				"python -m pip uninstall -y osprey JPype-py2"
			)
		}
	}

	val python3UninstallScripts by creating {
		group = "build"
		description = "Make uninstall scripts for python 3 distribution"
		doLast {
			writeScript(
				pythonBuildDir, "uninstall",
				"python -m pip uninstall -y osprey"
			)
		}
	}

	// insert some build steps before we build the python dist
	"python2DistZip" {
		if (python2 != null) {
			dependsOn(python2Wheel, makeDoc, python2InstallScripts, python2UninstallScripts)
		} else {
			enabled = false
		}
	}
	"python3DistZip" {
		if (python3 != null) {
			dependsOn(python3Wheel, makeDoc, python3InstallScripts, python3UninstallScripts)
		} else {
			enabled = false
		}
	}

	val testRunScript by creating {
		doLast {

			val classpath =
				(
					listOf("classes") +
					sourceSets["test"].runtimeClasspath
						.filter { it.extension == "jar" }
						.map { "lib/${it.name}" }
				)
				.joinToClasspath()

			writeScript(
				buildDir, "run",
				"java -cp \"$classpath\" $@"
			)
		}
	}

	"testDistZip" {
		dependsOn(testClasses, testRunScript)
	}

	val serviceTar by creating(Tar::class) {
		group = "distribution"
		description = "build the app server runtime for this version of the osprey service"
		dependsOn(jar)

		archiveBaseName.set(releaseNameService)
		archiveVersion.set(versionService)
		destinationDirectory.set(releasesDir.toFile())
		compression = Compression.BZIP2

		val dir = buildDir.resolve("service-$versionService")
		doFirst {
			dir.recreateFolder()

			// write the run script
			val libs = ArrayList<String>().apply {
				jar.outputs.files
					.forEach { add(it.name) }
				sourceSets["main"].runtimeClasspath
					.filter { it.extension == "jar" }
					.forEach { add(it.name) }
			}
			val classpath = libs.joinToClasspath { "lib/$it" }
			writeScript(
				dir, releaseNameService,
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
			from(sourceSets["main"].runtimeClasspath.filter { it.extension == "jar" })
		}
		into("progs") {
			from(projectDir.resolve("progs"))
		}

		// cleanup
		doLast {
			dir.deleteFolder()
		}
	}

	val serviceDockerTar by creating(Tar::class) {
		group = "distribution"
		description = "build the distribution package of the docker image for the osprey service"

		archiveBaseName.set(releaseNameServiceDocker)
		archiveVersion.set(versionService)
		destinationDirectory.set(releasesDir.toFile())

		// don't bother compressing this tar
		// the VAST MAJORITY of the space is taken up by the docker image, which is already compressed
		// we won't gain much more by compressing a few small text files
		//compression = Compression.BZIP2

		val imagePath = buildDir.resolve("docker/$releaseNameServiceDocker-$versionService.tar.bz2")
		val serviceDir = projectDir.resolve("src/main/docker/service")

		val dir = buildDir.resolve("service-docker")
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
			from("$serviceDir/$releaseNameService")
			from("$serviceDir/install.sh")
			from("$serviceDir/uninstall.sh")
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


	val archiveReleases by creating {
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
							releaseArchiveDir.resolve(release.fileName).toString(),
							SftpProgressLogger()
						)
					}
				} else {
					println("No new releases to upload")
				}
			}
		}
	}

	val downloadServiceReleases by creating {
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
							releaseArchiveDir.resolve(release.filename).toString(),
							releasesDir.resolve(release.filename).toString(),
							SftpProgressLogger()
						)
					}
				} else {
					println("No extra service releases to download")
				}
			}
		}
	}

	val compileShaders by creating {
		group = "build"
		doLast {

			val outDir = sourceSets["main"].resources.srcDirs.first()
				.resolve("$packagePath/molscope/shaders")

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

	this["build"].dependsOn(compileShaders)

	// add documentation to the gui distribution
	jpackageImage {
		doLast {
			val jp = project.runtime.jpackageData.get()
			val imageDir = when (os) {
				OperatingSystem.MAC_OS -> jp.imageOutputDir.resolve(jp.imageName + ".app").resolve("Contents")
				else -> jp.imageOutputDir.resolve(jp.imageName)
			}
			copy {
				from(
					projectDir.resolve("readme.md"),
					projectDir.resolve("LICENSE.txt")
				)
				into(imageDir)
			}
		}
	}

	// TODO: replace this with the nifty licenser plugin?
	//   https://github.com/Minecrell/licenser
	val updateLicenseHeaders by creating {
		group = "build"
		description = "updates license headers in all source files"
		doLast {
			updateLicenseHeaders()
		}
	}

	val javadocJsonFile = buildDir.resolve("doc/javadoc.json")

	val parseJavadoc by creating {
		group = "documentation"
		description = "export javadocs into a queryable format"
		dependsOn("testClasses")
		// NOTE: this task apparently won't re-run after a code recompile
		// my gradle-fu isn't good enough to figure out how to do that
		// in the meantime, just delete the build/doc/javadoc.json file to get this task to run again
		doLast {

			javadocJsonFile.parent.createFolder()
			javadocJsonFile.write { json ->

				javaexec {
					classpath = sourceSets.test.get().runtimeClasspath
					mainClass.set("build.JavadocTool")
					jvmArgs(
						*moduleArgs.toTypedArray()
					)
					args(
						packagePath.replace('/', '.'),
						sourceSets.main.get().java.sourceDirectories.first().absolutePath
					)
					standardOutput = json
				}
			}
		}
		outputs.file(javadocJsonFile)
	}

	// TODO: parse kotlin docs?

	val generatePythonDocs by creating {
		group = "documentation"
		dependsOn(parseJavadoc)
		inputs.files(javadocJsonFile)
		doLast {

			// generate the documentation into a hugo module in the build folder
			val modDir = buildDir.resolve("doc/code-python")
			modDir.recreateFolder()

			// init the hugo module
			exec {
				commandLine("hugo", "mod", "init", "code-python")
				workingDir = modDir.toFile()
			}
			val contentDir = modDir.resolve("content")
			contentDir.createFolder()

			// render python docs
			val modules = listOf(
				"osprey",
				"osprey.prep",
				"osprey.ccs",
				"osprey.slurm",
				"osprey.jvm"
			)
			for ((modulei, module) in modules.withIndex()) {
				pydocMarkdown(module, contentDir.resolve("$module.md"), weight = modulei + 1)
			}
		}
	}

	val generateJavaDocs by creating {
		group = "documentation"
		doLast {

			// generate the documentation into a hugo module in the build folder
			val modDir = buildDir.resolve("doc/code-java")
			modDir.recreateFolder()

			// init the hugo module
			exec {
				commandLine("hugo", "mod", "init", "code-java")
				workingDir = modDir.toFile()
			}
			val contentDir = modDir.resolve("content")
			contentDir.createFolder()

			// render java docs, see:
			// https://docs.oracle.com/en/java/javase/17/docs/specs/man/javadoc.html
			exec {
				commandLine(
					"javadoc",
					"-source", "$javaLangVersion",
					"-sourcepath", sourceSets.main.get().java.srcDirs.first().absolutePath,
					"-d", contentDir.toString(),
					"-subpackages", packagePath.replace('/', '.'),
					"-classpath", sourceSets.main.get().runtimeClasspath.joinToClasspath { it.absolutePath },
					*moduleArgs.toTypedArray(),
					"-Xdoclint:none"
				)
			}

			// hugo wants to use the index.html url,
			// so rename the index file generated by javadoc to something else
			contentDir.resolve("index.html").rename("start.html")

			// tweak the markdown files from the javadoc folder, otherwise hugo get confused
			Files.walk(contentDir).use { stream ->
				stream.asSequence()
					.filter { it.extension() == "md" }
					.forEach {
						val text = it.read()
						it.write {
							write(hugoFrontMatter(hidden = true))
							write(text)
						}
					}
			}
		}
	}

	val generateKotlinDocs by creating {
		group = "documentation"
		doLast {

			// generate the documentation into a hugo module in the build folder
			val modDir = buildDir.resolve("doc/code-kotlin")
			modDir.recreateFolder()

			// init the hugo module
			exec {
				commandLine("hugo", "mod", "init", "code-kotlin")
				workingDir = modDir.toFile()
			}
			val contentDir = modDir.resolve("content")
			contentDir.createFolder()

			// render Kotlin docs, see:
			// https://kotlin.github.io/dokka/1.5.30/user_guide/cli/usage/
			// https://discuss.kotlinlang.org/t/problems-running-dokka-cli-1-4-0-rc-jar-from-the-command-line/18855
			javaexec {
				mainClass.set("org.jetbrains.dokka.MainKt")
				setClasspath(sourceSets["test"].runtimeClasspath)
				args(
					"-moduleName", packagePath.replace('/', '.'),
					"-outputDir", contentDir.toString(),
					"-sourceSet", listOf(
						"-src", kotlin.sourceSets.main.get().kotlin.srcDirs.first().absolutePath,
						"-classpath", sourceSets.main.get().runtimeClasspath.joinToString(";") { it.absolutePath }
					).joinToString(" "),
					"-pluginsClasspath", sourceSets.test.get().runtimeClasspath
						// filter out the CLI jar from the plugin classpath... apparently dokka complains otherwise
						.filter { !it.nameWithoutExtension.startsWith("dokka-cli-") }
						.joinToString(";") { it.absolutePath }
					// NOTE: dokka always expects ; as the path separator, regardless of platform
				)
			}

			// hugo wants to use the index.html url,
			// so rename the index file generated by dokka to something else
			contentDir.resolve("index.html").rename("start.html")
		}
	}

	val generateCodeDocs by creating {
		group = "documentation"
		description = "Generate the Python, Java, and Kotlin code documentation for the current source tree"
		dependsOn(generatePythonDocs, generateJavaDocs, generateKotlinDocs)
	}

	fun checkHugoPrereqs() {

		// commands we'll need
		commandExistsOrThrow("hugo")
		commandExistsOrThrow("git")
		commandExistsOrThrow("go")

		// download the theme, if needed
		val themeDir = buildDir.resolve("doc/hugo-theme-learn")
		if (!themeDir.exists()) {

			exec {
				commandLine(
					"git", "clone",
					"--depth", "1",
					"--branch", "2.5.0",
					"https://github.com/matcornic/hugo-theme-learn",
					themeDir.toString()
				)
			}

			// pretend the theme is a go module
			exec {
				commandLine("go", "mod", "init", "local.tld/hugo-theme-learn")
				workingDir = themeDir.toFile()
			}
		}

		// make sure we got it
		if (!themeDir.exists()) {
			throw Error("Hugo theme is not available. The download must have failed somehow.")
		}
	}

	val buildWebsite by creating {
		group = "documentation"
		description = "Builds the Osprey documentation and download website"
		doLast {

			checkHugoPrereqs()

			val webDir = buildDir.resolve("website")
			webDir.recreateFolder()

			exec {
				commandLine(
					"hugo",
					"--destination", webDir.toString()
				)
				workingDir = docDir.toFile()
			}
		}
	}

	val hugoServer by creating {
		group = "documentation"
		description = "Start the Hugo development server. Useful for writing documentation"
		doLast {

			checkHugoPrereqs()

			val webDir = buildDir.resolve("website")

			exec {
				commandLine(
					"hugo",
					"server",
					"--destination", webDir.toString()
				)
				workingDir = docDir.toFile()
			}
		}
	}
}


/**
 * Hugo front matter, in TOML format, with learn theme extensions
 * https://gohugo.io/content-management/front-matter/
 * https://learn.netlify.app/en/cont/pages/#front-matter-configuration
 */
fun hugoFrontMatter(
	title: String? = null,
	menuTitle: String? = null,
	weight: Int = 4,
	disableToc: Boolean = false,
	hidden: Boolean = false
): String =
	"""
		|+++
		|${if (title != null) """title = "$title" """ else ""}
		|${if (menuTitle != null) """menuTitle = "$menuTitle" """ else ""}
		|weight = $weight
		|${if (disableToc) "disableToc = true" else ""} 
		|${if (hidden) "hidden = true" else ""} 
		|+++
		|
		|
	""".trimMargin()


fun pydocMarkdown(module: String, file: Path, title: String = module, weight: Int = 5) {

	// is pydoc-markdown installed
	commandExistsOrThrow("pydoc-markdown")

	val configPath = docDir.resolve("pydoc-markdown.yml")

	file.write { out ->

		// write the hugo front matter
		write(hugoFrontMatter(
			title,
			weight = weight,
			hidden = true
		))

		// flush buffers before pointing other external programs into this stream
		flush()

		// generate the markdown from the python module using pydoc-markdown
		// https://github.com/NiklasRosenstein/pydoc-markdown
		exec {
			commandLine(
				"pydoc-markdown",
				"--search-path", projectDir.resolve("src/main/python"),
				"--module", module,
				configPath.toString()
			)
			workingDir = projectDir.toFile()
			standardOutput = out
			environment["PYTHONPATH"] = listOf(
				System.getenv("PYTHONPATH"),
				projectDir.resolve("src/test/python")
			).joinToClasspath()
		}
	}
}

fun nvcc(exec: Exec, kernelName: String, maxRegisters: Int? = null, profile: Boolean = false) {

	val args = mutableListOf("nvcc")

	if (profile) {
		// if profiling, compile for one arch with profiling/debug info
		// NOTE: change this to your GPU's arch
		args.addAll(listOf("-cubin", "-gencode=arch=compute_61,code=sm_61", "-lineinfo", "--ptxas-options=-v"))
	} else {
		// otherwise, compile for all archs
		// see Maxwell compatibility guide:
		// http://docs.nvidia.com/cuda/maxwell-compatibility-guide/index.html#building-maxwell-compatible-apps-using-cuda-6-0
		args.addAll(listOf("-fatbin",
			"-gencode=arch=compute_20,code=sm_20",
			"-gencode=arch=compute_30,code=sm_30",
			"-gencode=arch=compute_35,code=sm_35",
			"-gencode=arch=compute_50,code=sm_50",
			"-gencode=arch=compute_52,code=sm_52",
			"-gencode=arch=compute_60,code=sm_60",
			"-gencode=arch=compute_61,code=sm_61",
			"-gencode=arch=compute_62,code=sm_62",
			"-gencode=arch=compute_62,code=compute_62"
		))
	}

	if (maxRegisters != null) {
		args.addAll(listOf("-maxrregcount", "$maxRegisters"))
	}

	args.addAll(listOf("$kernelName.cu", "-o", "$kernelName.bin"))

	exec.workingDir = file("src/main/resources/gpuKernels/cuda")
	exec.commandLine(args)
}

fun String.writeToFile(path: Path, newline: String) {
	BufferedWriter(OutputStreamWriter(Files.newOutputStream(path), StandardCharsets.UTF_8.newEncoder())).use { out ->
		for (line in split("\n")) {
			out.append(line)
			out.append(newline)
		}
	}
}

fun writeScript(dir: Path, filename: String, cmd: String) {
	when (os) {
		OperatingSystem.MAC_OS,
		OperatingSystem.LINUX -> {

			val file = dir.resolve("$filename.sh")

			"""
				|#! /bin/sh
				|$cmd
			""".trimMargin()
			.writeToFile(file, "\n")

			// set the shell script executable
			Files.setPosixFilePermissions(file, Files.getPosixFilePermissions(file).apply {
				add(PosixFilePermission.OWNER_EXECUTE)
			})
		}
		OperatingSystem.WINDOWS -> {

			val file = dir.resolve("$filename.bat")

			"""
				|@echo off
				|$cmd
			""".trimMargin()
			.writeToFile(file, "\r\n")
		}
		else -> throw Error("unrecognized operating system: $os")
	}
}


enum class HeaderResult {
	Updated,
	Ignored
}


// TODO: replace this custom license header code with the licenser plugin
// https://github.com/Minecrell/licenser

fun updateLicenseHeaders() {

	val header = """
		|This file is part of OSPREY 3.0
		|
		|OSPREY Protein Redesign Software Version 3.0
		|Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
		|
		|OSPREY is free software: you can redistribute it and/or modify
		|it under the terms of the GNU General Public License version 2
		|as published by the Free Software Foundation.
		|
		|You should have received a copy of the GNU General Public License
		|along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
		|
		|OSPREY relies on grants for its development, and since visibility
		|in the scientific literature is essential for our success, we
		|ask that users of OSPREY cite our papers. See the CITING_OSPREY
		|document in this distribution for more information.
		|
		|Contact Info:
		|   Bruce Donald
		|   Duke University
		|   Department of Computer Science
		|   Levine Science Research Center (LSRC)
		|   Durham
		|   NC 27708-0129
		|   USA
		|   e-mail: www.cs.duke.edu/brd/
		|
		|<signature of Bruce Donald>, Mar 1, 2018
		|Bruce Donald, Professor of Computer Science
		""".trimMargin()
		.lines()

	val deleteTheseAutoHeaders = listOf(
		""" |/*
			| * To change this template, choose Tools | Templates
			| * and open the template in the editor.
			| */
		""".trimMargin(),
		""" |/*
			| * To change this license header, choose License Headers in Project Properties.
			| * To change this template file, choose Tools | Templates
			| * and open the template in the editor.
			| */
		""".trimMargin()
	)

	fun applyCHeader(lines: MutableList<String>): HeaderResult {

		// extract the existing license header, if any
		var readMode = 0
		var existingHeader = lines
			.takeWhile {
				val line = it.trim()
				when (readMode) {
					0 -> {
						if (line.startsWith("/*")) {
							readMode = 1
							return@takeWhile true
						}
					}
					1 -> {
						if (line.startsWith("**")) {
							return@takeWhile true
						} else if (line.startsWith("*/")) {
							readMode = 2
							return@takeWhile true
						}
					}
				}
				return@takeWhile false
			}
			.map { it.substring(2).trim() }
		if (existingHeader.size >= 3) {
			for (i in 0 until existingHeader.size) {
				lines.removeAt(0)
			}
			existingHeader = existingHeader.subList(1, existingHeader.size - 1)
		}

		// if it matches the desired header, then we're done
		if (existingHeader == header) {
			return HeaderResult.Ignored
		}

		// trim blank lines
		while (lines.firstOrNull()?.isBlank() == true) {
			lines.removeAt(0)
		}

		// add the new header
		lines.add(0, "")
		lines.add(0, "*/")
		for (i in 0 until header.size) {
			lines.add(0, "** " + header[header.size - i - 1])
		}
		lines.add(0, "/*")

		return HeaderResult.Updated
	}

	fun applyPythonHeader(lines: MutableList<String>): HeaderResult {

		// extract the existing license header, if any
		val existingHeader = lines
			.takeWhile { it.startsWith("##") }
			.map { it.substring(2).trim() }
		for (i in 0 until existingHeader.size) {
			lines.removeAt(0)
		}

		// if it matches the desired header, then we're done
		if (existingHeader == header) {
			return HeaderResult.Ignored
		}

		// trim blank lines
		while (lines.firstOrNull()?.isBlank() == true) {
			lines.removeAt(0)
		}

		// add the new header
		lines.add(0, "")
		for (i in 0 until header.size) {
			lines.add(0, "## " + header[header.size - i - 1])
		}

		return HeaderResult.Updated
	}

	fun applyHeader(path: Path, applier: (MutableList<String>) -> HeaderResult) {

		var text = Files.readAllBytes(path).toString(StandardCharsets.UTF_8)

		// remove any headers automatically added by IDEs or other tools
		var removedAutoHeaders = false
		for (autoHeader in deleteTheseAutoHeaders) {
			val newtext = text.replace(autoHeader, "")
			if (newtext != text) {
				removedAutoHeaders = true
				text = newtext
			}
		}

		val lines = text.lines().toMutableList()

		// trim blank lines from the top
		while (lines.firstOrNull()?.isBlank() == true) {
			lines.removeAt(0)
		}

		// keep one blank line on the bottom
		// NOTE: a trailing newline creates a blank line at the end of the file,
		// so it's sufficient to remove all blank entries in the lines list
		while (lines.lastOrNull()?.isBlank() == true) {
			lines.removeAt(lines.size - 1)
		}

		// anything left?
		if (lines.isEmpty()) {
			return
		}

		if (removedAutoHeaders || applier(lines) == HeaderResult.Updated) {
			Files.write(path, lines)
			println("updated: $path")
		}
	}

	fun applyHeaders(dirname: String, filter: (String) -> Boolean, applier: (MutableList<String>) -> HeaderResult) {

		// for each matched file in the folder (and subfolders)
		val dir = projectDir.resolve(dirname)
		Files.walk(dir)
			.filter { filter(it.fileName.toString()) }
			.forEach { applyHeader(it, applier) }
	}

	// apply header to java files
	for (dirname in listOf("src", "test")) {
		applyHeaders(
			dirname,
			filter = { it.endsWith(".java") },
			applier = ::applyCHeader
		)
	}

	// apply header to this file
	applyHeader(projectDir.resolve("build.gradle.kts"), ::applyCHeader)

	// apply header to kernel files
	for (dirname in listOf("src/main/resources/gpuKernels")) {
		applyHeaders(
			dirname,
			filter = { it.endsWith(".cu") || it.endsWith(".cl") },
			applier = ::applyCHeader
		)
	}

	// apply header to python files
	for (dirname in listOf("python")) {
		applyHeaders(
			dirname,
			filter = { it.endsWith(".py") },
			applier = ::applyPythonHeader
		)
	}

	// NOTE: don't apply the header to the python example scripts.
	// there's no need to scare osprey users with legalese in the tutorials
}


fun Project.propertyOrNull(name: String): Any? =
	try {
		property(name)
	} catch (ex: groovy.lang.MissingPropertyException) {
		null
	}


/** configure and open an SFTP connection over SSH */
fun <T> sftp(block: ChannelSftp.() -> T): T {

	// get user auth info
	val sshDir = Paths.get(System.getProperty("user.home")).resolve(".ssh")
	val user = project.propertyOrNull("dlab.user") as String?
		?: throw Error("no user configured. set `dlab.user = your-user` in gradle.properties")
	val keypriv = (project.propertyOrNull("dlab.key.private") as String?)
		?.let { Paths.get(it) }
		?: sshDir.resolve("id_rsa")
	val keypub = (project.propertyOrNull("dlab.key.public") as String?)
		?.let { Paths.get(it) }
		?: Paths.get("$keypriv.pub")
	val keytype = (project.propertyOrNull("dlab.key.type") as String?)
	val knownHosts = (project.propertyOrNull("dlab.knownHosts") as String?)
		?.let { Paths.get(it) }
		?: sshDir.resolve("known_hosts")

	// configure host info
	val host = "login.cs.duke.edu"
	val port = 22

	// API docs: http://www.jcraft.com/jsch/
	val jsch = JSch()
	jsch.addIdentity(keypriv.toString(), keypub.toString(), null)
	jsch.setKnownHosts(knownHosts.toString())

	// if the key type is known, set it explicitly
	// since older SSH server versions don't support key type negotiation,
	// and sometimes SSH clients hit the auth failure threshold before finding the right algorithm via brute force, see:
	// https://github.com/mwiede/jsch/issues/45#issuecomment-839727926
	if (keytype != null) {
		JSch.setConfig("PubkeyAcceptedAlgorithms", keytype)
	}

	// capture the JSch log
	val log = StringBuilder()
	JSch.setLogger(object : JschLogger {
		override fun isEnabled(level: Int) = true
		override fun log(level: Int, message: String?) {
			val levelstr = when (level) {
				JschLogger.DEBUG -> "DEBUG"
				JschLogger.INFO -> "INFO"
				JschLogger.WARN -> "WARN"
				JschLogger.ERROR -> "ERROR"
				JschLogger.FATAL -> "FATAL"
				else -> "???"
			}
			log.append("\t$levelstr: $message\n")
		}
	})

	// open the SSH connection
	val session = jsch.getSession(user, host, port)
	session.setConfig("PreferredAuthentications", "publickey")
	try {
		session.connect()
	} catch (t: Throwable) {
		System.err.println("""
			|Error connecting to SSH server. Troubleshooting tips:
			|   Make sure your username is correct: $user  (the username should not be quoted)
			|   Check that the private key is correct: $keypriv (exists? ${Files.exists(keypriv)})
			|      If your private key is not at the default location, set `dlab.key.private = /path` in gradle.properties.
			|   Check that the public key is correct: $keypub (exists? ${Files.exists(keypub)})
			|      If your public key is not at the default location, set `dlab.key.public = /path` in gradle.properties.
			|   Make sure the SSH client can find your key type before the auth failure limit (sometimes as low as 2).
			|      Check the SSH log (below) for details.
			|      If the correct key type isn't found before the auth failure limit, try setting `dlab.key.type = your-key=type` in gradle.properties.
			|      For older SSH keys, the correct type is often `ssh-rsa`.
			|   Check that the known hosts file is correct: $knownHosts (exists? ${Files.exists(knownHosts)})
			|      If your known_hosts file is not at the default location, set `dlab.knownHosts = /path` in gradle.properties.
			|   Make sure the SSH host is in your known_hosts file. Try connecting to $host via `ssh` in a terminal first.
			|   You'll need a Duke CS account to connect to the Duke CS SSH server.
			|SSH log:
			|$log
		""".trimMargin())
		throw t
	}
	try {
		val channel = session.openChannel("sftp") as ChannelSftp
		try {
			channel.connect()

			// at long last, we're connected. Do The Thing!
			return channel.block()

		} finally {
			channel.disconnect()
		}
	} finally {
		if (session.isConnected) {
			session.disconnect()
		}
	}
}

class SftpProgressLogger : SftpProgressMonitor {

	private var bytes: Long? = null
	private var progress: Long = 0
	private var startNs: Long? = null
	private var lastNs: Long? = null

	override fun init(op: Int, src: String?, dest: String?, max: Long) {
		println("Copying $max bytes from $src to $dest ...")
		bytes = max
		startNs = System.nanoTime()
		lastNs = startNs
	}

	override fun count(count: Long): Boolean {

		progress += count

		// don't log more than once per second
		val nowNs = System.nanoTime()
		val lastNs = lastNs
		val elapsedNs = if (lastNs != null) {
			nowNs - lastNs
		} else {
			0
		}
		if (elapsedNs >= 1_000_000_000) {

			this.lastNs = nowNs

			val bytes = bytes
			if (bytes != null) {
				val pct = 100f*progress.toFloat()/bytes.toFloat()
				println("   copied $progress bytes ${"%.1f".format(pct)} %")
			} else {
				println("   copied $progress bytes")
			}
		}

		return true // true to continue, false to cancel
	}

	override fun end() {
		val nowNs = System.nanoTime()
		val startNs = startNs
		if (startNs != null) {
			println("   done in ${"%.1f".format((nowNs - startNs).toFloat()/1e9f)} s")
		} else {
			println("   done")
		}
	}
}


// some conveniences for files and paths

fun Path.exists(): Boolean =
	Files.exists(this)

fun Path.deleteFolder() =
	toFile().deleteRecursively()

fun Path.createFolder(): Path {
	Files.createDirectories(this)
	return this
}

fun Path.recreateFolder() {
	deleteFolder()
	createFolder()
}

fun Path.copyTo(to: Path) =
	Files.copy(this, to)

fun Path.copyFolderTo(to: Path) =
	toFile().copyRecursively(to.toFile(), overwrite = true)

fun Path.listFiles(): Sequence<Path> =
	Files.list(this).asSequence()
		.filter { Files.isRegularFile(it) }

fun Path.read(): String =
	Files.readString(this)

fun Path.write(block: java.io.Writer.(java.io.OutputStream) -> Unit) {
	toFile().outputStream().use { out ->
		out.writer().use { writer ->
			writer.block(out)
		}
	}
}

fun Path.rename(dst: String) {
	if ('/' in dst) {
		throw IllegalArgumentException("invalid rename, can't have /")
	}
	Files.move(this, parent.resolve(dst))
}

fun Path.extension(): String? =
	fileName.toString().split('.')
		.takeIf { it.size > 1 }
		?.last()

fun Path.deleteFile() =
	Files.delete(this)


/**
 * Returns true iff the command exists
 */
fun commandExists(name: String): Boolean {

	// check the cache first
	commandExistence[name]?.let { return it }

	val exists = when (os) {

		OperatingSystem.WINDOWS -> {
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

fun commandExistsOrThrow(name: String) {
	if (!commandExists(name)) {
		throw Error("command not available: $name")
	}
}


fun <T> Iterable<T>.joinToClasspath(transform: ((T) -> CharSequence)? = null): String {
	if (transform == null) {
		return joinToString(File.pathSeparator)
	}

	return this.map(transform).joinToString(File.pathSeparator)
}
