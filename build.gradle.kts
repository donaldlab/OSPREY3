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
import java.net.URL
import java.nio.charset.StandardCharsets
import java.nio.file.attribute.PosixFilePermission
import java.nio.file.StandardCopyOption
import java.nio.file.Files
import java.nio.file.Path
import org.gradle.internal.os.OperatingSystem
import org.jetbrains.kotlin.gradle.tasks.KotlinCompile



plugins {
	`java-library`
	kotlin("jvm") version "1.3.60"
	kotlin("plugin.serialization") version "1.3.61"
	application
	idea
	id("org.openjfx.javafxplugin") version("0.0.7")
	id("org.beryx.runtime") version "1.8.0"
}

javafx {
	version = "11"
	modules("javafx.controls")
}

val projectDir = project.projectDir.toPath().toAbsolutePath()
val pythonSrcDir = projectDir.resolve("src/main/python")
val pythonBuildDir = projectDir.resolve("build/python")
val pythonWheelDir = pythonBuildDir.resolve("wheel")
val pythonWheelhouseDir = pythonSrcDir.resolve("wheelhouse")
val docSrcDir = pythonSrcDir.resolve("doc")
val docBuildDir = pythonBuildDir.resolve("doc")


group = "edu.duke.cs"
version = Files.readAllLines(projectDir.resolve("src/main/resources/config/version"))[0]

val versionService = "0.3"

repositories {
	jcenter()
}

java {
	sourceCompatibility = JavaVersion.VERSION_14
	targetCompatibility = JavaVersion.VERSION_14
}

idea {
	module {
		// use the same output folders as gradle, so the pythonDevelop task works correctly
		outputDir = sourceSets["main"].output.classesDirs.singleFile
		testOutputDir = sourceSets["test"].output.classesDirs.singleFile
		inheritOutputDirs = false
	}
}

application {
	mainClassName = "edu.duke.cs.osprey.control.Main"
}

dependencies {

	// kotlin runtime
	implementation(kotlin("stdlib-jdk8"))

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
	implementation("ch.obermuhlner:big-math:2.0.1")
	implementation("org.joml:joml:1.9.19")
	implementation("org.tukaani:xz:1.8")
	implementation("com.hazelcast:hazelcast:4.0")
	implementation("net.java.dev.jna:jna:5.5.0")
	implementation("com.google.guava:guava:29.0-jre")

	// TODO: get rid of "api" dependencies after merging everything into the majestic monolith

	// libs used by the GUI
	api("org.apache.commons:commons-lang3:3.4")
	api("commons-io:commons-io:2.5")
	api("org.tomlj:tomlj:1.0.0")
	api(files("lib/kdtree.jar")) // no authoritative source on the internet

	// used by the gui
	api("com.cuchazinteractive:kludge:0.1")

	// used by the service
	val ktorVersion = "1.3.0"
	api("io.ktor:ktor-server-netty:$ktorVersion")
	implementation("io.ktor:ktor-serialization:$ktorVersion")

	// for JCuda, gradle tries (and fails) download the natives jars automatically,
	// so turn off transitive dependencies. we'll deal with natives manually
	implementation("org.jcuda:jcuda:0.8.0") {
		isTransitive = false
	}

	// build systems never seem to like downloading from arbitrary URLs for some reason...
	// so make a helper func to do it for us
	fun url(urlString: String): ConfigurableFileCollection {

		// parse the URL
		val url = URL(urlString)
		val filename = urlString.split("/").last()
		val libDir = projectDir.resolve("lib")
		val filePath = libDir.resolve(filename)

		// create a gradle task to download the file
		val downloadTask by tasks.creating {
			if (!Files.exists(filePath)) {
				Files.createDirectories(libDir)
				url.openStream().use {
					Files.copy(it, filePath, StandardCopyOption.REPLACE_EXISTING)
				}
			}
		}

		// return a files collection that depends on the task
		return files(filePath) {
			builtBy(downloadTask)
		}
	}

	// TPIE-Java isn't in the maven/jcenter repos yet, download directly from Github
	// use `api` instead of `implementation` here, since the compiler sometimes complains about these types not being available in downstream projects
	api(url("https://github.com/donaldlab/TPIE-Java/releases/download/v1.1/edu.duke.cs.tpie-1.1.jar"))

	// native libs for GPU stuff
	listOf("natives-linux-amd64", "natives-macosx-universal", "natives-windows-amd64").forEach {
		runtimeOnly("org.jogamp.gluegen:gluegen-rt:2.3.2:$it")
		runtimeOnly("org.jogamp.jocl:jocl:2.3.2:$it")
	}
	listOf("linux-x86_64", "apple-x86_64", "windows-x86_64").forEach {
		runtimeOnly("org.jcuda:jcuda-natives:0.8.0:$it")
	}
}

// add the module dependencies directly to the javac args
// I don't think gradle has a good way to handle this yet?
val moduleArgs = listOf(
	"--add-modules=jdk.incubator.foreign"
)
// TODO: add the module args for distributions too?

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
		val classpathStr = classpath.joinToString(File.pathSeparatorChar.toString())
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

		jvmTarget = "1.8"

		// enable experimental features so we can use the fancy ktor stuff
		freeCompilerArgs += "-Xuse-experimental=kotlin.Experimental"
		freeCompilerArgs += "-XXLanguage:+InlineClasses"
	}
}

tasks.withType<Test> {
	// the default 512m is too little memory to run test designs
	maxHeapSize = "2g"
	useJUnitPlatform()
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
		"jdk.zipfs" // needed to provide jar:// file system
	)
}


val os = OperatingSystem.current()

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

// get the OS name
val osname: String by lazy {
	when (os) {
		OperatingSystem.LINUX -> "linux"
		OperatingSystem.WINDOWS -> "win"
		OperatingSystem.MAC_OS -> "osx"
		else -> throw Error("unrecognized operating system: $os")
	}
}


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

	// tell gradle to write down the version numbers where apps can read them
	processResources {

		// always update the build properties
		outputs.upToDateWhen { false }

		// write down the osprey service version
		from(sourceSets["main"].resources.srcDirs) {
			include("edu/cs/duke/osprey/molscope/build.properties")
			expand(
				"dev" to isDev
			)
		}

		// write down the osprey service version
		from(sourceSets["main"].resources.srcDirs) {
			include("edu/cs/duke/osprey/service/build.properties")
			expand(
				"version" to versionService
			)
		}
	}

	val appendBuildNumber by creating {
		dependsOn("processResources")
		doLast {

			// read the version number
			val versionFile = projectDir.resolve("build/resources/main/config/version").toFile()
			var version = versionFile.readText().trim()

			// append the CI build ID, if available

			version += if (rootProject.hasProperty("AZURE_BUILD_ID")) {
				val versionId = rootProject.property("AZURE_BUILD_ID")
				".$versionId\n"
			} else {
				"-dev\n"
			}
			versionFile.writeText(version)
		}
	}

	val jar = "jar" {
		dependsOn(appendBuildNumber)
	}.get()

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
		dependsOn("runtime")

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
				.joinToString(":")

			writeScript(
				buildDir.toPath(), "run",
				"java -cp \"$classpath\" $@"
			)
		}
	}

	"testDistZip" {
		dependsOn(testClasses, testRunScript)
	}

	val compileShaders by creating {
		group = "build"
		doLast {

			val outDir = sourceSets["main"].resources.srcDirs.first()
				.resolve("edu/duke/cs/osprey/molscope/shaders")

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

	// TODO: replace this with the nifty licenser plugin?
	//   https://github.com/Minecrell/licenser
	val updateLicenseHeaders by creating {
		group = "build"
		description = "updates license headers in all source files"
		doLast {
			updateLicenseHeaders()
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
